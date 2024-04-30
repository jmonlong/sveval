#include <Rcpp.h>
using namespace Rcpp;

//' Prepare clusters from pairwise overlap information. Assumes that only one overlap
//' relationship is present for each pair, i.e. only "SV1 overlaps SV2", no
//' "SV2 overlaps SV1". Also no self-overlap. This is done by keeping only
//' "queryHits<subjectHits", for example.
//' @title Prepare clusters from pairwise overlap information
//' @param queryHits the IDs of the first SVs
//' @param subjectHits the IDs of the second SVs
//' @param total total number of SVs
//' @return list with two vectors: 'cl' with cluster ID, 'clique' with boolean about clique configuration
//' @author Jean Monlong
//' @keywords internal
// [[Rcpp::export]]
List prepare_cluster_cpp(IntegerVector queryHits, IntegerVector subjectHits, IntegerVector totalNbSVs){
    std::vector<int> cls(totalNbSVs[0], 0);
    std::vector<bool> cliques(totalNbSVs[0], true);
    std::map<int, int> edge_per_cl;
    // loop over the edge list and update cluster assignments
    for(unsigned int eidx=0; eidx<queryHits.size(); eidx++){
        if((cls[queryHits[eidx]-1] == 0) & (cls[subjectHits[eidx]-1] == 0)){
            // none of the SVs have a cluster assigned yet
            // assign them to a new cluster named as the smaller index
            int cl_idx = subjectHits[eidx];
            if(queryHits[eidx] < subjectHits[eidx]){
                cl_idx = queryHits[eidx];
            }
            cls[queryHits[eidx]-1] = cl_idx;
            cls[subjectHits[eidx]-1] = cl_idx;
            // init the edge count to 1
            edge_per_cl[cl_idx] = 1;
        } else if((cls[queryHits[eidx]-1] != 0) & (cls[subjectHits[eidx]-1] != 0)){
            // both have a cluster already assigned
            // if they're not the same clusters, merge them
            if(cls[queryHits[eidx]-1] != cls[subjectHits[eidx]-1]){
                // assign them to the cluster with the smaller id
                int cl_idx = cls[subjectHits[eidx]-1];
                int ocl_idx = cls[queryHits[eidx]-1];
                if(ocl_idx < cl_idx){
                    cl_idx = cls[queryHits[eidx]-1];
                    ocl_idx = cls[subjectHits[eidx]-1];
                }
                // change all the other cluster's assignments
                for(unsigned int oidx=0; oidx<cls.size(); oidx++){
                    if(cls[oidx] == ocl_idx){
                        cls[oidx] = cl_idx;
                    }
                }
                // update the edge counts too
                edge_per_cl[cl_idx] += edge_per_cl[ocl_idx] + 1;
            } else {
                // they are already assigned to the same cluster
                // just update the edge tally
                edge_per_cl[cls[queryHits[eidx]-1]]++;
            }
        } else if(cls[queryHits[eidx]-1] != 0){
            // queryHits SV has a cluster assigned, assign the same to subjectHits
            cls[subjectHits[eidx]-1] = cls[queryHits[eidx]-1];
            // update the edge counts
            edge_per_cl[cls[queryHits[eidx]-1]]++;
        } else if(cls[subjectHits[eidx]-1] != 0){
            // queryHits SV has a cluster assigned, assign the same to subjectHits
            cls[queryHits[eidx]-1] = cls[subjectHits[eidx]-1];
            // update the edge counts
            edge_per_cl[cls[subjectHits[eidx]-1]]++;
        }
    }
    // compute the size of each cluster
    std::map<int, int> cl_size;
    for(unsigned int idx=0; idx<cls.size(); idx++){
        if(cl_size.find(cls[idx]) != cl_size.end()){
            cl_size[cls[idx]]++;
        } else {
            cl_size[cls[idx]] = 1;
        }
    }
    // fill in the SVs without clusters and guess which one is a clique
    for(unsigned int idx=0; idx<cls.size(); idx++){
        if(cls[idx] == 0){
            cls[idx] = idx + 1;
        } else {
            if(2 * edge_per_cl[cls[idx]] < (cl_size[cls[idx]] * (cl_size[cls[idx]] - 1))){
                // cliques should n(n-1)/2 edges (n being the size of the cluster)
                // if not flag as not a clique
                cliques[idx] = false;
            }
        }
    }
    return List::create(_["cl"] = cls, _["clique"] = cliques);
}
