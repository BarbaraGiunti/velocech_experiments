#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <chrono>
#include <utility> // Required for std::pair
#include <algorithm>
#include <cmath>
#include <limits>
#include <phat_helper_functions.h>
#include <thread>
#include <array>
#include "inputs.h"
#ifdef _OPENMP
#include <parallel/algorithm>
#endif
using namespace phat;

class edge;

std::vector<std::vector<double>> pointCloud;
std::vector<std::vector<uint16_t>> nborsSortedByDist;    // for each point p, the other points sorted by distance to p
std::vector<int32_t> edgeIndices;                          // edgeIndices[pairIndex(a, b)] = position of edge ab in edges

// std::vector<double> edgeDistances;
int dim;
std::vector<double> coords; // flat copy of pointCloud: point a is coords[a*dim .. a*dim+dim-1]
inline double sqDist(index a, index b){
    const double* x = &coords[a*dim];
    const double* y = &coords[b*dim];
    double s = 0;
    for(int l=0; l<dim; l++){
        double d = x[l]-y[l];
        s += d*d;
    }
    return s;
}
int32_t noOfPoints;
index noOfPairs;
// position of the pair {a, b} (a != b) in the list of all pairs i<j, ordered by i and then j
inline index pairIndex(index a, index b){
    if(a > b) std::swap(a, b);
    return a*noOfPoints - a*(a+1)/2 + (b - a - 1);
}
bit_tree_column tempCol;
index H1pivotcounter = 0;
double lowerThresh = 0;

std::vector<edge> edges;
std::vector<column> cols;                  // columns of edges with edgeInfo1 == -2
#define COLM(e) (cols[(e).edgeInfo2])      // for those edges, edgeInfo2 is the slot in cols

class edge
{
    protected:
    public:
        double distance;
        double deathTime;
        int32_t edgeInfo1;
        int32_t edgeInfo2;
        std::pair<uint16_t, uint16_t> points;
        
        // should only be called if edgeInfo1 >= -1, which is the same as i being crowded
        void simplifyNonImp(){
            bool madeAChange = true;
            while(madeAChange){
                madeAChange = false;
                if(edgeInfo1 >= 0 && edges[edgeInfo1].deathTime <= std::max(distance, lowerThresh)){
                    edge& edgeA = edges[edgeInfo1];
                    if(edgeA.edgeInfo1 == -1){
                        edgeInfo1 = edgeA.edgeInfo2;
                        madeAChange = true;
                    }
                    else if(edgeA.edgeInfo1 == -2){
                        if(COLM(edgeA).size() == 0){
                            edgeInfo1 = -1;
                            madeAChange = true;
                        }
                        if(COLM(edgeA).size() == 1){
                            edgeInfo1 = COLM(edgeA)[0];
                            madeAChange = true;
                        }
                    }
                    else if(edgeA.edgeInfo2 == -1){
                        edgeInfo1 = edgeA.edgeInfo1;
                        madeAChange = true;
                    }
                    else if(edgeInfo2 == -1 && edgeA.edgeInfo1 >= -1){
                        edgeInfo1 = edgeA.edgeInfo1;
                        edgeInfo2 = edgeA.edgeInfo2;
                        madeAChange = true;
                    }
                }
                if(edgeInfo2 >= 0 && edges[edgeInfo2].deathTime <= std::max(distance, lowerThresh)){
                    edge& edgeB = edges[edgeInfo2];
                    if(edgeB.edgeInfo1 == -1){
                        edgeInfo2 = edgeB.edgeInfo2;
                        madeAChange = true;
                    }
                    else if(edgeB.edgeInfo1 == -2){
                        if(COLM(edgeB).size() == 0){
                            edgeInfo2 = -1;
                            madeAChange = true;
                        }
                        if(COLM(edgeB).size() == 1){
                            edgeInfo2 = COLM(edgeB)[0];
                            madeAChange = true;
                        }
                    }
                    else if(edgeB.edgeInfo2 == -1){
                        edgeInfo2 = edgeB.edgeInfo1;
                        madeAChange = true;
                    }
                    else if(edgeInfo1 == -1 && edgeB.edgeInfo1 >= -1){
                        edgeInfo1 = edgeB.edgeInfo1;
                        edgeInfo2 = edgeB.edgeInfo2;
                        madeAChange = true;
                    }
                }
            }
        }
        //returns true iff edge i is important
        bool computemu(index pIndex, index qIndex){
            for (index m : nborsSortedByDist[pIndex]){
                if (m == qIndex) {
                    edgeInfo1 = -3;
                    edgeInfo2 = 0;   // must not be -1 for uncrowded edges (see simplifyNonImp)
                    return true;
                }
                if(sqDist(pIndex, m) + sqDist(m, qIndex) < distance) {
                    edgeInfo1 = edgeIndices[pairIndex(m, pIndex)];
                    edgeInfo2 = edgeIndices[pairIndex(m, qIndex)];
                    return false;
                }
            }
            return true;
        }

        // should only be called for a == -2, which means i is important and its column has been initialized
        // assumes columns are stored without the max el.
        void simplifyImp(){
            if(COLM(*this).empty()){
                return;
            }
            const edge& edgeM = edges[COLM(*this).back()];
            // b is dead, or is important but not initialized, or
            // column of b is to the right of that of i, so we cannot add it.
            if(edgeM.edgeInfo1 == -3 || edgeM.deathTime > std::max(deathTime, lowerThresh)){
                return;
            }
            COLM(*this).pop_back();
            if(edgeM.edgeInfo1 == -1 && edgeM.edgeInfo2 == -1){
                return;
            }
            tempCol.add_col(COLM(*this));
            if(edgeM.edgeInfo1 == -2){
                tempCol.add_col(COLM(edgeM));
            }
            else{
                if(edgeM.edgeInfo1 >= 0){
                    tempCol.add_index(edgeM.edgeInfo1);
                }
                if(edgeM.edgeInfo2 >= 0){
                    tempCol.add_index(edgeM.edgeInfo2);
                }
            }
            tempCol.get_col_and_clear(COLM(*this));
        }
};

static_assert(sizeof(edge) == 32, "edge should be 32 bytes");

int findIndex(index p, const std::vector<uint16_t>& nbrs, double target){
    auto it = std::upper_bound(nbrs.begin(), nbrs.end(), target,
        [p](double value, uint16_t r) { return value < sqDist(p, r); });
    if (it == nbrs.begin()) {
        return -1;  // No element <= query
    }
    return static_cast<int>(it - nbrs.begin()) - 1;
}


// returns true if the bitcol is 0 in the end
// If doNotTouchOtherColumns, then return true if the bitcol is 0 in the end
// If !doNotTouchOtherColumns, then return true if some column has gotten a new deathTime
// (or whatever the precise statement is)
bool reduceImp(bit_tree_column* bitCol, double& deathT, bool doNotTouchOtherColumns, int inst){
    bool returnVal = false;
    while(!bitCol->is_empty()){
        index m = bitCol->get_max_index();
        // m not important
        if(edges[m].edgeInfo1 >= -1 && edges[m].deathTime <= deathT){
            bitCol->add_index(m);
            edges[m].simplifyNonImp();
            if(edges[m].edgeInfo1 >= 0){
                bitCol->add_index(edges[m].edgeInfo1);
            }
            if(edges[m].edgeInfo2 >= 0){
                bitCol->add_index(edges[m].edgeInfo2);
            }
            continue;
        }
        // add column of m normally
        if(edges[m].edgeInfo1 == -2 && edges[m].deathTime <= deathT){
            bitCol->add_index(m);
            edges[m].simplifyImp();
            bitCol->add_col(COLM(edges[m]));
            continue;
        }
        if(doNotTouchOtherColumns){
            return false;
        }
        bitCol->add_index(m);
        // found new pivot
        if(edges[m].edgeInfo1 == -3){
            edges[m].deathTime = deathT;
            column col;
            (*bitCol).get_col_and_clear(col);
            if(col.size() == 0){
                edges[m].edgeInfo1 = -1;
                edges[m].edgeInfo2 = -1;
            }
            else if(col.size() == 1){
                edges[m].edgeInfo1 = col[0];
                edges[m].edgeInfo2 = -1;
            }
            else if(col.size() == 2){
                edges[m].edgeInfo1 = col[0];
                edges[m].edgeInfo2 = col[1];
            }
            else{
                edges[m].edgeInfo1 = -2;
                edges[m].edgeInfo2 = (int32_t)cols.size();
                cols.push_back(std::move(col));
            }
            return true;
        }
        // swap with column that has later deathTime
        if(edges[m].deathTime > deathT){
            returnVal = true;
            const double temp = edges[m].deathTime;
            edges[m].deathTime = deathT;
            deathT = temp;
            column tempcol;
            bitCol->get_col(tempcol);
            const bool hadCol = (edges[m].edgeInfo1 == -2);
            if(hadCol){
                bitCol->add_col(COLM(edges[m]));
            }
            else{
                if(edges[m].edgeInfo1 >= 0){
                    bitCol->add_index(edges[m].edgeInfo1);
                }
                if(edges[m].edgeInfo2 >= 0){
                    bitCol->add_index(edges[m].edgeInfo2);
                }
            }
            if(tempcol.size() > 2){
                if(!hadCol){                       // needs a new slot
                    edges[m].edgeInfo2 = (int32_t)cols.size();
                    cols.emplace_back();
                }
                edges[m].edgeInfo1 = -2;
                std::swap(COLM(edges[m]), tempcol);
            }
            else{
                if(hadCol){
                    column().swap(COLM(edges[m])); // free the old column before edgeInfo2 is overwritten
                }
                if(tempcol.size() == 0){
                    edges[m].edgeInfo1 = -1;
                    edges[m].edgeInfo2 = -1;
                }
                else if(tempcol.size() == 1){
                    edges[m].edgeInfo1 = tempcol[0];
                    edges[m].edgeInfo2 = -1;
                }
                else{
                    edges[m].edgeInfo1 = tempcol[0];
                    edges[m].edgeInfo2 = tempcol[1];
                }
            }
        }
    }
    if(!doNotTouchOtherColumns){
        return returnVal;
    }
    return true;
}

double threshForZeroingOut(index i){
    double thresh = -INFINITY;
    tempCol.clear();
    tempCol.add_index(i);
    while(!tempCol.is_empty()){
        index m = tempCol.get_max_index();
        // m crowded
        if(edges[m].edgeInfo1 >= -1){
            tempCol.add_index(m);
            if(edges[m].edgeInfo1 >= 0){
                tempCol.add_index(edges[m].edgeInfo1);
            }
            if(edges[m].edgeInfo2 >= 0){
                tempCol.add_index(edges[m].edgeInfo2);
            }
            thresh = std::max(thresh, edges[m].deathTime);
            continue;
        }
        // add column of m normally
        if(edges[m].edgeInfo1 == -2){
            tempCol.add_index(m);
            tempCol.add_col(COLM(edges[m]));
            thresh = std::max(thresh, edges[m].deathTime);
            continue;
        }
        tempCol.clear();
        return INFINITY;
    }
    tempCol.clear();
    return thresh;
}

//One long vector with an element for each pair i<j of point indices, sorted by sq distance.
void allSortedDistances() {
    edges.resize(noOfPairs);
    #pragma omp parallel for schedule(dynamic, 16)
    for (index i=0; i < noOfPoints; i++){
        // position of pair (i, i+1) in the list of all pairs i<j
        index ijCounter = i*noOfPoints - i*(i+1)/2;
        for (index j=i+1; j < noOfPoints; j++){
            edge& e = edges[ijCounter];
            e.distance = sqDist(i, j);
            e.deathTime = -1;
            e.points = {(uint16_t)i, (uint16_t)j};
            ijCounter ++;
        }
    }
    auto byDist = [](const edge& lhs, const edge& rhs) { return lhs.distance < rhs.distance; };
    #ifdef _OPENMP
        // balanced quicksort sorts in place; the default parallel sort needs a full extra copy of edges
        __gnu_parallel::sort(edges.begin(), edges.end(), byDist, __gnu_parallel::balanced_quicksort_tag());
    #else
        std::sort(edges.begin(), edges.end(), byDist);
    #endif
    edgeIndices.resize(noOfPairs);
    #pragma omp parallel for
    for (index k=0; k < noOfPairs; k++){
        edgeIndices[pairIndex(edges[k].points.first, edges[k].points.second)] = (int32_t)k;
    }
}

// Claude's new version
// For each point p, the list of all other points sorted by distance to p.
std::vector<std::vector<uint16_t>> sortedDistancesFromEachPoint() {
    std::vector<std::vector<uint16_t>> distFromEachPoint(noOfPoints);
    for(auto& v : distFromEachPoint) v.reserve(noOfPoints - 1);
    for(const edge& e : edges){
        const auto [a, b] = e.points;
        distFromEachPoint[a].push_back(b);
        distFromEachPoint[b].push_back(a);
    }
    return distFromEachPoint;
}

std::vector<std::pair<double, double>> H1barcode(){
    std::vector<std::pair<index, index>> columnsH0(noOfPoints, std::pair{-1, -1});
    std::vector<edge*> H1PivotEdges(0);
    for(index k=0; k < noOfPairs; k++){
        edge& edg = edges[k];
        const auto [p, q] = edg.points;
        if(!edg.computemu(p,q)){
            // edge k is crowded
            continue;
        }
        int pivot = q;
        int nonpivot = p;
        while(true){
            //If the pivot is already a pivot for an earlier column
            if(columnsH0[pivot].second != -1){
                int otherNonpivot = columnsH0[pivot].first;
                if(nonpivot == otherNonpivot){
                    //the column goes to zero; this means edge #i contributes to H1
                    H1pivotcounter++;
                    H1PivotEdges.push_back(&edg);
                    break;
                }
                pivot = std::max(nonpivot, otherNonpivot);
                nonpivot = std::min(nonpivot, otherNonpivot);
            }
            else{
                //haven't seen this pivot before, so store the reduced column and move on
                columnsH0[pivot] = {nonpivot, pivot};
                //the edge is dead (so can be ignored when computing H1)
                edg.edgeInfo1 = -1;
                edg.edgeInfo2 = -1;
                break;
            }
        }
    }
    /////////////////////////////////////
    ////         Computing H1        ////
    /////////////////////////////////////
    column zeroColumn = {};
    std::vector<std::pair<double, double>> barcode(H1pivotcounter);
    bit_tree_column bitColumnSingleEdge, bitColumnpqr;
    bitColumnSingleEdge.init(noOfPairs);
    bitColumnpqr.init(noOfPairs);
    index i = -1;
    for(edge& edg : edges){
        i++;
        double pq = edg.distance;
        bitColumnSingleEdge.add_index(i);
        if(reduceImp(&bitColumnSingleEdge, pq, true, 1)){
            edg.edgeInfo1 = -1;
            edg.edgeInfo2 = -1;
            continue;
        }
        const index pIndex = edg.points.first;
        const index qIndex = edg.points.second;
        column columnSingleEdge = {};
        bitColumnSingleEdge.get_col_and_clear(columnSingleEdge);
        double threshold = threshForZeroingOut(i);
        double diamThreshold = 0;
        ////////////////////////////////////////////
        ////  Go through the cofaces of edge i  ////
        ////////////////////////////////////////////
        const std::vector<uint16_t>& nborsP = nborsSortedByDist[pIndex];
        for(index j=0; j < (index)nborsP.size(); j++){
            const index rIndex = nborsP[j];
            if(rIndex == qIndex){
                j = std::max(j, (index)findIndex(pIndex, nborsP, diamThreshold-pq));
                continue;
            }
            //////////////////////////////////////////////////////////////////////////
            ////  A bunch of geometric checks that may allow skipping iterations  ////
            //////////////////////////////////////////////////////////////////////////
            const double pr = sqDist(pIndex, rIndex);
            const double qr = sqDist(qIndex, rIndex);
            // If the triangle is obtuse, it's crowded, so we can skip it
            if(pq + pr <= qr || pq + qr <= pr){
                continue;
            }
            // this is the square of the diameter of the circumcircle of pqr
            double pqrDiameter = 2*pq*pr*qr/(pq*pr+pr*qr+qr*pq-(pq*pq+pr*pr+qr*qr)/2);
            if(pr + qr <= pq){
                diamThreshold = std::max(diamThreshold, pqrDiameter);
                if(threshold < diamThreshold - 0.00001){
                    break;
                }
                continue;
            }
            if(pqrDiameter < diamThreshold - 0.00001){
                continue;
            }
            if(pr > threshold){
                break;
            }
            if(pqrDiameter > threshold){
                continue;
            }
            ////////////////////////////////////////////////////////////
            ////  Done with geometry, the rest is column reduction  ////
            ////////////////////////////////////////////////////////////
            // bitColumnpqr is empty
            index prIndex = edgeIndices[pairIndex(pIndex, rIndex)];
            index qrIndex = edgeIndices[pairIndex(qIndex, rIndex)];
            bitColumnpqr.add_index(prIndex);
            bitColumnpqr.add_index(qrIndex);
            bitColumnpqr.add_col(columnSingleEdge);
            if(reduceImp(&bitColumnpqr, pqrDiameter, false, 0)){
                threshold = threshForZeroingOut(i);
                if(threshold < diamThreshold - 0.00001){
                    break;
                }
            }
        }
    }
    for(index i=0; i < H1pivotcounter; i++){
        barcode[i] = {(*H1PivotEdges[i]).distance, (*H1PivotEdges[i]).deathTime};
    }
    return barcode;
}

int main(int argc, char** argv) {
    std::string input_filename, output_filename;
    parse_command_line(argc, argv, input_filename, output_filename);
    if(!load_ascii(input_filename, pointCloud, dim)){
        return 1;
    }
    if(pointCloud.size() < 2 || pointCloud.size() > 65536){
        std::cerr << "Error: need between 2 and 65536 points, got " << pointCloud.size() << "\n";
        return 1;
    }
    noOfPoints = (int32_t)pointCloud.size();
    noOfPairs = (index)noOfPoints*(noOfPoints-1)/2;
    tempCol.init(noOfPairs);

    // Timing goes to std::cerr, so that standard output contains only the barcode.
    std::cerr << noOfPoints << " points in dimension " << dim << "\n";
    coords.resize((size_t)noOfPoints*dim);
    for(int a=0; a<noOfPoints; a++) for(int l=0; l<dim; l++) coords[(size_t)a*dim+l] = pointCloud[a][l];
    auto start = std::chrono::steady_clock::now();
    allSortedDistances();
    nborsSortedByDist = sortedDistancesFromEachPoint();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> firstTime = end - start;
    std::cerr << "Time for distances: " << firstTime.count() << " seconds\n";
    auto start2 = std::chrono::steady_clock::now();
    std::vector<std::pair<double, double>> barcode = H1barcode();
    auto end2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> secondTime = end2 - start2;
    std::cerr << "Time after computing distances: " << secondTime.count() << " seconds\n";

    // Bars are (birth^2, death^2); take square roots for the actual barcode.
    if(!write_barcode(output_filename, barcode)){
        return 1;
    }
    return 0;
}