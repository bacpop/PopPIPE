import os, sys
import numpy as np
import subprocess
import tempfile
import dendropy

def update_distance_matrices(refList, distMat, queryList = None, query_ref_distMat = None,
                             query_query_distMat = None):
    """Convert distances from long form (1 matrix with n_comparisons rows and 2 columns)
    to a square form (2 NxN matrices), with merging of query distances if necessary.

    Args:
        refList (list)
            List of references
        distMat (numpy.array)
            Two column long form list of core and accessory distances
            for pairwise comparisons between reference db sequences
        queryList (list)
            List of queries
        query_ref_distMat (numpy.array)
            Two column long form list of core and accessory distances
            for pairwise comparisons between queries and reference db
            sequences
        query_query_distMat (numpy.array)
            Two column long form list of core and accessory distances
            for pairwise comparisons between query sequences

    Returns:
        seqLabels (list)
            Combined list of reference and query sequences
        coreMat (numpy.array)
            NxN array of core distances for N sequences
        accMat (numpy.array)
            NxN array of accessory distances for N sequences
    """
    seqLabels = refList
    if queryList is not None:
        seqLabels = seqLabels + queryList

    coreMat = np.zeros((len(seqLabels), len(seqLabels)))
    accMat = np.zeros((len(seqLabels), len(seqLabels)))

    # Fill in symmetric matrices for core and accessory distances
    i = 0
    j = 1

    # ref v ref (used for --create-db)
    for row in distMat:
        coreMat[i, j] = row[0]
        coreMat[j, i] = coreMat[i, j]
        accMat[i, j] = row[1]
        accMat[j, i] = accMat[i, j]

        if j == len(refList) - 1:
            i += 1
            j = i + 1
        else:
            j += 1

    # if query vs refdb (--assign-query), also include these comparisons
    if queryList is not None:

        # query v query - symmetric
        i = len(refList)
        j = len(refList)+1
        for row in query_query_distMat:
            coreMat[i, j] = row[0]
            coreMat[j, i] = coreMat[i, j]
            accMat[i, j] = row[1]
            accMat[j, i] = accMat[i, j]
            if j == (len(refList) + len(queryList) - 1):
                i += 1
                j = i + 1
            else:
                j += 1

        # ref v query - asymmetric
        i = len(refList)
        j = 0
        for row in query_ref_distMat:
            coreMat[i, j] = row[0]
            coreMat[j, i] = coreMat[i, j]
            accMat[i, j] = row[1]
            accMat[j, i] = accMat[i, j]
            if j == (len(refList) - 1):
                i += 1
                j = 0
            else:
                j += 1

    # return outputs
    return seqLabels, coreMat, accMat

with open(snakemake.input['pkl'], 'rb') as pickle_file:
    rlist, qlist, self = pickle.load(pickle_file)
complete_distMat = np.load(snakemake.input['npy'])
combined_seq, core_distMat, acc_distMat = update_distance_matrices(rlist, complete_distMat)

tmp_core = tempfile.NamedTemporaryFile(delete=False)
phylip_name = tmp_core.name
tree_filename = snakemake.output[0]

with open(phylip_name, 'w') as pFile:
    pFile.write(str(len(rlist))+"\n")
    for coreDist, ref in zip(core_distMat, rlist):
        pFile.write(ref)
        pFile.write(' '+' '.join(map(str, coreDist)))
        pFile.write("\n")

# construct tree
rapidnj_cmd = "rapidnj " + phylip_name + " -i pd -o t -x " + tree_filename + ".raw"
try:
    # run command
    subprocess.run(rapidnj_cmd, shell=True, check=True)

    # remove quotation marks for microreact
    with open(tree_filename + ".raw", 'r') as f, open(tree_filename, 'w') as fo:
        for line in f:
            fo.write(line.replace("'", ''))
    os.remove(tree_filename+".raw")

# record errors
except subprocess.CalledProcessError as e:
    sys.stderr.write("Could not run command " + rapidnj_cmd + "; returned code: " + str(e.returncode) + "\n")
    sys.exit(1)

os.remove(phylip_name)