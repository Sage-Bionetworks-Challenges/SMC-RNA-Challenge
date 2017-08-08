#include "SetPair.h"

int numIntersect(vector<int> & s1, vector<int> & s2)
{

  vector<int> v(s1.size()+s2.size());
  vector<int>::iterator it;

  sort (s1.begin(),s1.end());  
  sort (s2.begin(),s2.end());   

  it=set_intersection (s1.begin(), s1.end(), s2.begin(), s2.end(), v.begin());
                                               
  v.resize(it-v.begin()); 
  
  return v.size();
}

vector<int> getUnion(vector<int> & s1, vector<int> & s2)
{
  vector<int> v(s1.size()+s2.size());
  vector<int>::iterator it;

  sort (s1.begin(),s1.end());
  sort (s2.begin(),s2.end());

  it=set_union (s1.begin(), s1.end(), s2.begin(), s2.end(), v.begin());

  v.resize(it-v.begin());

  return v;
}


int SetPair::merge(vector<set_pair_t> & spv)
{
    if(spv.size()==0)
        return 0;
    bool changed=true;
    while(changed==true)
    {
        changed=false;
        vector<int> keep(spv.size(),1);
        for(int i=0;i<spv.size()-1;i++)
        {
            for(int j=i+1;j<spv.size();j++)
            {
                if(numIntersect(spv[i].ids1, spv[j].ids1)>0 && numIntersect(spv[i].ids2, spv[j].ids2)>0)
                {
                    keep[i]=0;
                    spv[j].ids1=getUnion(spv[i].ids1,spv[j].ids1);
                    spv[j].ids2=getUnion(spv[i].ids2,spv[j].ids2);
                    changed=true;
                    break;
                }
            }
        }
        vector<set_pair_t> spv2;
        for(int i=0;i<spv.size();i++)
        {
            if(keep[i]==1)
                spv2.push_back(spv[i]);
        }
        spv=spv2;
    }
    return 0;
}

//use after merge
int SetPair::numIntersectSet(vector<set_pair_t> & res,vector<set_pair_t> & truth, evaluate_t & et)
{
    vector<int> found(truth.size(),0);
    vector<int> correct(res.size(),0);

    for(int i=0;i<res.size();i++)
    {
        for(int j=0;j<truth.size();j++)
        {
            int inter1=numIntersect(res[i].ids1, truth[j].ids1);
            int inter2=numIntersect(res[i].ids2, truth[j].ids2);
            if(inter1>0 && inter2>0) {
                found[j]=1;
                correct[i]=1;
                et.gene_true_positives.push_back(truth[j]);
            }
        }
    }

    int num=0;
    for(int i=0;i<found.size();i++) {
      if (found[i] == 1) {
        num += 1;
      } else {
        et.gene_false_negatives.push_back(truth[i]);
      }
    }

    for(int i=0;i<correct.size();i++) {
      if (correct[i] == 0) {
        et.gene_false_positives.push_back(res[i]);
      }
    }

    return num;
}




