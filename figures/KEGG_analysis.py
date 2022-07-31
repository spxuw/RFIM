import sklearn.feature_selection as skf
import numpy as np
import gseapy
import mygene
import pdb
#import subprocess

# adopted from the amim-test-suit at: https://github.com/dbblumenthal/amim-test-suite/

def compute_neg_log_gsea_p_value(pathways, result_genes):
    """Computes gene set enrichment score for result genes w.r.t. phenotype-related pathways.
    Parameters
    ----------
    pathways : list of str
        Names of phenotype-related pathways.
    result_genes : list of str
        Set of genes computed by network enrichment algorithm.
    Returns
    -------
    mean neg_log_gsea_p_value : float
        Negative log-transformed p-value of gene set enrichment analysis.
    """
    if len(result_genes) == 0:
        return 0.0
    try:
        mg = mygene.MyGeneInfo()

        out = mg.querymany(result_genes, scopes= 'entrezgene', fields='symbol', species='human', verbose=False)
        #gene_names = []
        #for line in out:
        #    try:
        #        gene_names.append(line["symbol"])
        #    except KeyError:
        #        pass
        res = gseapy.enrichr(gene_list=result_genes, description='pathway', gene_sets='KEGG_2016', cutoff=0.05,
                             outdir='../temp/enrichment', no_plot=True)
        full_results = res.results
        terms = list(full_results.Term)
        terms = [x.split(' ')[-1] for x in terms]
        p_values = []
        for i in range(len(terms)):
            if terms[i] in pathways:
                p_values.append(-np.log10(full_results['Adjusted P-value'][i]))
        if len(p_values) > 0:
            return np.mean(p_values)
        else:
            return 0
    except:
        return -1

# KEGG of true networks
phenotype = ['Asthma','Breast_cancer','Lung_cancer','Colorectal_cancer','Gastric_cancer',
              'Prostate_cancer','CVD_new','Diabetes_new']
random_type = ['RDPN','rewired','scale_free','shuffled','uniform']
random_type = ['shuffled']

pathways = {
	"Asthma": "hsa05310",
	"Breast_cancer": "hsa05224",
	"Lung_cancer": ["hsa05222","hsa05223"],
	"Colorectal_cancer": "hsa05210",
	"Gastric_cancer": "hsa05226",
	"Prostate_cancer": "hsa05215",
	"CVD_new": ["hsa05410","hsa05412","hsa05414","hsa05416"],
	"Diabetes_new": ["hsa04940","hsa04930"]
}

for k in range(8):
	gene_true = np.loadtxt('../../results/RFIM/network-'+str(phenotype[k])+'_String_f_0'+'.txt-spin1.txt',dtype=str)
	gene_true = gene_true.tolist()
	p_true = compute_neg_log_gsea_p_value(pathways[phenotype[k]],gene_true)
	for i in range(1):
		p_array = []
		p_array.append(p_true)
		for j in range(10):
			print(j)
			gene_rand = np.loadtxt('../../results/RFIM/network-'+str(phenotype[k])+'_String_'+str(random_type[i])+'_'+str(j)+'.txt-spin1.txt',dtype=str)
			gene_rand = gene_rand.tolist()
			p_rand = compute_neg_log_gsea_p_value(pathways[phenotype[k]],gene_rand)
			p_array.append(p_rand)
		np.savetxt('../../results/GSEA_meaningfulness/'+str(phenotype[k])+'_'+random_type[i]+'_String.txt',p_array)

for k in range(8):
    gene_true = np.loadtxt('../../results/RFIM/network-'+str(phenotype[k])+'_iRefIndex_f_0'+'.txt-spin1.txt',dtype=str)
    gene_true = gene_true.tolist()
    p_true = compute_neg_log_gsea_p_value(pathways[phenotype[k]],gene_true)
    for i in range(1):
        p_array = []
        p_array.append(p_true)
        for j in range(10):
            print(j)
            gene_rand = np.loadtxt('../../results/RFIM/network-'+str(phenotype[k])+'_iRefIndex_'+str(random_type[i])+'_'+str(j)+'.txt-spin1.txt',dtype=str)
            gene_rand = gene_rand.tolist()
            p_rand = compute_neg_log_gsea_p_value(pathways[phenotype[k]],gene_rand)
            p_array.append(p_rand)
        np.savetxt('../../results/GSEA_meaningfulness/'+str(phenotype[k])+'_'+random_type[i]+'_iRefIndex.txt',p_array)

# other methods
methods = ['DIAMOnD','ModuleDiscoverer','DOMINO','ROBUST','Hierarchical HotNet']
for i in range(5):
	for k in range(8):
		p_array=[]
		gene_true = np.loadtxt('../../results/module_txt/'+str(methods[i])+'_'+str(phenotype[k])+'_String.txt',dtype=str)
		gene_true = gene_true.tolist()
		p_true = compute_neg_log_gsea_p_value(pathways[phenotype[k]],gene_true)
		p_array.append(p_true)
		np.savetxt('../../results/GSEA_meaningfulness/'+str(methods[i])+'_'+str(phenotype[k])+'_String.txt',p_array)


for i in range(5):
	for k in range(8):
		p_array=[]
		gene_true = np.loadtxt('../../results/module_txt/'+str(methods[i])+'_'+str(phenotype[k])+'_iRefIndex.txt',dtype=str)
		gene_true = gene_true.tolist()
		p_true = compute_neg_log_gsea_p_value(pathways[phenotype[k]],gene_true)
		p_array.append(p_true)
		np.savetxt('../../results/GSEA_meaningfulness/'+str(methods[i])+'_'+str(phenotype[k])+'_iRefIndex.txt',p_array)
