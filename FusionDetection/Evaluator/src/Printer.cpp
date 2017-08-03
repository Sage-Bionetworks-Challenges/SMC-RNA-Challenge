#include "Printer.h"

int Printer::print(vector<evaluate_t> & evals, char * outputfile,  string gene_file, int base_resolution, int max_diff, string pseudo_counts, string version, Gene g)
{
    ofstream ofs;
    ofs.open (outputfile, std::ofstream::out);

    ofs<<"#fusionToolEvaluator version "+version<<endl;
    ofs<<"#--gene-annotation-file="+gene_file+"; --base-resolution="+my_db2string(base_resolution)+"; --max-diff="+my_db2string(max_diff)+"; --pseudo-counts="+pseudo_counts<<endl;
    ofs<<"#Evaluation_Name\tNum_Res_Trans\tNum_Truth_Trans\tSensitivity_Trans\tPrecision_Trans\tF1_Trans\tNum_Res_Gene\tNum_Truth_Gene\tSensitivity_Gene\tPrecision_Gene\tF1_Gene"<<endl;
    for(int i=0;i<evals.size();i++)
    {
        evaluate_t et=evals[i];
        ofs<<et.name<<"\t";
        ofs<<et.num_res_trans<<"\t";
        ofs<<et.num_truth_trans<<"\t";
        ofs<<et.sensitivity_t<<"\t";
        ofs<<et.precision_t<<"\t";
        ofs<<et.f_t<<"\t";
        ofs<<et.num_res_gene<<"\t";
        ofs<<et.num_truth_gene<<"\t";
        ofs<<et.sensitivity_g<<"\t";
        ofs<<et.precision_g<<"\t";
        ofs<<et.f_g<<"\n";
	ofs<<"\n";
	for (int i = 0; i < et.false_negatives.size(); i++) {
            ofs<<"FN: " + et.false_negatives[i].name<<"\n";
        }
	ofs<<"\n";
	for (int i = 0; i < et.false_positives.size(); i++) {
            ofs<<"FP: " + et.false_positives[i].name<<"\n";
        }
	ofs<<"\n";
	for (int i = 0; i < et.true_positives.size(); i++) {
            ofs<<"TP: " + et.true_positives[i].name<<"\n";
        }
	ofs<<"\n";

	for (int i = 0; i < et.gene_true_positives.size(); i++) {
    set_pair_t sp = et.gene_true_positives[i];
    for (int j = 0; j < sp.ids1.size(); j++) {
      gene_t *x = g.getGene(sp.ids1[j]);
      for (int k = 0; k < sp.ids2.size(); k++) {
        gene_t *y = g.getGene(sp.ids2[k]);
        ofs<<"G_TP: " + x->name2 + "-" + y->name2<<"\n";
      }
    }
	  ofs<<"\n";
  }
	ofs<<"\n";

	for (int i = 0; i < et.gene_false_positives.size(); i++) {
    set_pair_t sp = et.gene_false_positives[i];
    for (int j = 0; j < sp.ids1.size(); j++) {
      gene_t *x = g.getGene(sp.ids1[j]);
      for (int k = 0; k < sp.ids2.size(); k++) {
        gene_t *y = g.getGene(sp.ids2[k]);
        ofs<<"G_FP: " + x->name2 + "-" + y->name2<<"\n";
      }
    }
	  ofs<<"\n";
  }

	ofs<<"\n";
	for (int i = 0; i < et.gene_false_negatives.size(); i++) {
    set_pair_t sp = et.gene_false_negatives[i];
    for (int j = 0; j < sp.ids1.size(); j++) {
      gene_t *x = g.getGene(sp.ids1[j]);
      for (int k = 0; k < sp.ids2.size(); k++) {
        gene_t *y = g.getGene(sp.ids2[k]);
        ofs<<"G_FN: " + x->name2 + "-" + y->name2<<"\n";
      }
    }
	  ofs<<"\n";
  }
	ofs<<"\n";

  /*
	ofs<<"resSP\n";
            ofs<<et.resSP[i].ids1+"-"+et.resSP[i].ids2<<"\n";
        }

	ofs<<"\n";
	ofs<<"truthSP\n";
	for (int i = 0; i < et.truthSP.size(); i++) {
            ofs<<et.truthSP[i].ids1+"-"+et.truthSP[i].ids2<<"\n";
        }
    */
    }
    ofs.close();

    return 0;
}
