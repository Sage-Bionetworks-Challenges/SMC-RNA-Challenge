#include "Printer.h"

int Printer::print(vector<evaluate_t> & evals, string outputfile,  string gene_file, int base_resolution, int max_diff, string pseudo_counts, string version, Gene g)
{
    ofstream ofs;
    ofs.open ((char *)outputfile.c_str(), std::ofstream::out);

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

        ofstream ofs_fuse;
        ofs_fuse.open ((char *)(outputfile + std::to_string(i) + ".fusions").c_str(), std::ofstream::out);

        // Dump transcript fusion results (TP, FP, FN)
        for (int i = 0; i < et.false_negatives.size(); i++)
        {
            std::string n = et.false_negatives[i].name;
            n = n.replace(n.find("-"), 1, "\t");
            ///std::replace(n.begin(), n.end(), "-", "\t");
            ofs_fuse<<"T_FN\t" + n <<"\n";
        }

        for (int i = 0; i < et.false_positives.size(); i++)
        {
            std::string n = et.false_positives[i].name;
            n = n.replace(n.find("-"), 1, "\t");
            //std::replace(n.begin(), n.end(), '-', '\t');
            ofs_fuse<<"T_FP\t" + n <<"\n";
        }

        for (int i = 0; i < et.true_positives.size(); i++)
        {
            std::string n = et.true_positives[i].name;
            n = n.replace(n.find("-"), 1, "\t");
            //std::replace(n.begin(), n.end(), '-', '\t');
            ofs_fuse<<"T_TP\t" + n <<"\n";
        }

        for (int i = 0; i < et.gene_true_positives.size(); i++)
        {
          set_pair_t sp = et.gene_true_positives[i];
          for (int j = 0; j < sp.ids1.size(); j++)
          {
            gene_t *x = g.getGene(sp.ids1[j]);
            for (int k = 0; k < sp.ids2.size(); k++)
            {
              gene_t *y = g.getGene(sp.ids2[k]);
              ofs_fuse<<"G_TP\t" + x->name2 + "\t" + y->name2<<"\n";
            }
          }
        }

        for (int i = 0; i < et.gene_false_positives.size(); i++)
        {
          set_pair_t sp = et.gene_false_positives[i];
          for (int j = 0; j < sp.ids1.size(); j++)
          {
            gene_t *x = g.getGene(sp.ids1[j]);
            for (int k = 0; k < sp.ids2.size(); k++)
            {
              gene_t *y = g.getGene(sp.ids2[k]);
              ofs_fuse<<"G_FP\t" + x->name2 + "\t" + y->name2<<"\n";
            }
          }
        }

        for (int i = 0; i < et.gene_false_negatives.size(); i++)
        {
          set_pair_t sp = et.gene_false_negatives[i];
          for (int j = 0; j < sp.ids1.size(); j++)
          {
            gene_t *x = g.getGene(sp.ids1[j]);
            for (int k = 0; k < sp.ids2.size(); k++)
            {
              gene_t *y = g.getGene(sp.ids2[k]);
              ofs_fuse<<"G_FN\t" + x->name2 + "\t" + y->name2<<"\n";
            }
          }
        }

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
        ofs_fuse.close();
    }
    ofs.close();

    return 0;
}
