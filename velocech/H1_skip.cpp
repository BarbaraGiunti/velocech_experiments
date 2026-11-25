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
#include <cstring>
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
// For dim >= 8: four independent partial sums, so consecutive additions do not have to wait for each other.
// Optimization for very high-dimensional point clouds, where distance computation is a bottleneck.
// Kept out of line so that the low-dimensional code path below is compiled exactly as before.
[[gnu::noinline]] double sqDistLong(const double* x, const double* y){
    double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
    int l = 0;
    for(; l + 4 <= dim; l += 4){
        const double d0 = x[l]-y[l], d1 = x[l+1]-y[l+1], d2 = x[l+2]-y[l+2], d3 = x[l+3]-y[l+3];
        s0 += d0*d0; s1 += d1*d1; s2 += d2*d2; s3 += d3*d3;
    }
    for(; l < dim; l++){
        const double d = x[l]-y[l];
        s0 += d*d;
    }
    return (s0 + s1) + (s2 + s3);
}
inline double sqDist(index a, index b){
    const double* x = &coords[a*dim];
    const double* y = &coords[b*dim];
    if(dim >= 8) return sqDistLong(x, y);
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
// Death times are only stored for edges that can die: the uncrowded edges that create an H1 class.
// All other edges (crowded, or killed in H0) have death time -1, as before.
std::vector<double> deathTimes;
std::vector<int32_t> slotEdge;             // for each slot, the position of its edge in edges
std::vector<column> cols;                  // columns of edges with edgeInfo1 == -2
#define COLM(e) (cols[(e).edgeInfo2])      // for those edges, edgeInfo2 is the slot in cols

class edge
{
    protected:
    public:
        int32_t edgeInfo1;                     // during sorting, edgeInfo1 and edgeInfo2 hold the sort key
        int32_t edgeInfo2;
        int32_t slot;                          // index in deathTimes, or -1 for edges whose death time is always -1
        std::pair<uint16_t, uint16_t> points;
        double distance() const { return sqDist(points.first, points.second); }
        double deathTime() const { return slot < 0 ? -1.0 : deathTimes[slot]; }
        void setDeath(double d){ deathTimes[slot] = d; }
        
        // should only be called if edgeInfo1 >= -1, which is the same as i being crowded
        void simplifyCrowded(){
            double myDist = -1;   // max(distance, lowerThresh), computed the first time it is needed
            auto validNow = [&](double t){
                if(t < 0) return true;
                if(myDist < 0) myDist = std::max(distance(), lowerThresh);
                return t <= myDist;
            };
            bool madeAChange = true;
            while(madeAChange){
                madeAChange = false;
                if(edgeInfo1 >= 0 && validNow(edges[edgeInfo1].deathTime())){
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
                if(edgeInfo2 >= 0 && validNow(edges[edgeInfo2].deathTime())){
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
            const double distance = sqDist(pIndex, qIndex);
            for (index m : nborsSortedByDist[pIndex]){
                if (m == qIndex) {
                    edgeInfo1 = -3;
                    edgeInfo2 = 0;   // must not be -1 for uncrowded edges (see simplifyCrowded)
                    return true;
                }
                if(sqDist(pIndex, m) + sqDist(m, qIndex) < distance) {
                    edgeInfo1 = edgeIndices[pairIndex(m, pIndex)];
                    edgeInfo2 = edgeIndices[pairIndex(m, qIndex)];
                    return false;
                }
            }
            edgeInfo1 = -3;
            edgeInfo2 = 0;
            return true;
        }

        // should only be called for a == -2, which means i is important and its column has been initialized
        // assumes columns are stored without the max el.
        void simplifyUnc(){
            if(COLM(*this).empty()){
                return;
            }
            const edge& edgeM = edges[COLM(*this).back()];
            // b is dead, or is important but not initialized, or
            // column of b is to the right of that of i, so we cannot add it.
            if(edgeM.edgeInfo1 == -3 || edgeM.deathTime() > std::max(deathTime(), lowerThresh)){
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

// Used by exhaust().
bit_tree_column exhCol;

// col represents an edge at all times >= d (a column about to be stored for an edge that dies at time d).
// Rewrites col so that none of its entries can be expanded at time d: every entry that is crowded, zero,
// or an uncrowded edge that died at time <= d is replaced by its own representation, recursively.
// Each remaining entry is an uncrowded edge that is alive or dies later than d.
// Such columns are added many times later (always at times >= d), and without this the same
// expansions would be repeated every time. The result is in increasing order, like get_col.
void exhaust(column& col, double d){
    if(col.empty()) return;
    exhCol.add_col(col);
    column kept;
    while(!exhCol.is_empty()){
        const index m = exhCol.get_max_index();
        exhCol.add_index(m);   // remove m
        edge& em = edges[m];
        if(em.edgeInfo1 >= -1 && em.deathTime() <= d){
            em.simplifyCrowded();
            if(em.edgeInfo1 >= 0) exhCol.add_index(em.edgeInfo1);
            if(em.edgeInfo2 >= 0) exhCol.add_index(em.edgeInfo2);
        }
        else if(em.edgeInfo1 == -2 && em.deathTime() <= d){
            em.simplifyUnc();
            exhCol.add_col(COLM(em));
        }
        else{
            kept.push_back(m);   // the largest remaining entry is final: later expansions only add smaller indices
        }
    }
    std::reverse(kept.begin(), kept.end());
    col.swap(kept);
}

static_assert(sizeof(edge) == 16, "edge should be 16 bytes");

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
bool reduceUnc(bit_tree_column* bitCol, double& deathT, bool doNotTouchOtherColumns, int inst){
    bool returnVal = false;
    while(!bitCol->is_empty()){
        index m = bitCol->get_max_index();
        // m not important
        if(edges[m].edgeInfo1 >= -1 && edges[m].deathTime() <= deathT){
            bitCol->add_index(m);
            edges[m].simplifyCrowded();
            if(edges[m].edgeInfo1 >= 0){
                bitCol->add_index(edges[m].edgeInfo1);
            }
            if(edges[m].edgeInfo2 >= 0){
                bitCol->add_index(edges[m].edgeInfo2);
            }
            continue;
        }
        // add column of m normally
        if(edges[m].edgeInfo1 == -2 && edges[m].deathTime() <= deathT){
            bitCol->add_index(m);
            edges[m].simplifyUnc();
            bitCol->add_col(COLM(edges[m]));
            continue;
        }
        if(doNotTouchOtherColumns){
            return false;
        }
        bitCol->add_index(m);
        // found new pivot
        if(edges[m].edgeInfo1 == -3){
            edges[m].setDeath(deathT);
            column col;
            (*bitCol).get_col_and_clear(col);
            exhaust(col, deathT);
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
        if(edges[m].deathTime() > deathT){
            returnVal = true;
            const double temp = edges[m].deathTime();
            edges[m].setDeath(deathT);
            deathT = temp;
            column tempcol;
            bitCol->get_col(tempcol);
            exhaust(tempcol, edges[m].deathTime());   // m's new representation is valid from its new death time
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
            thresh = std::max(thresh, edges[m].deathTime());
            continue;
        }
        // add column of m normally
        if(edges[m].edgeInfo1 == -2){
            tempCol.add_index(m);
            tempCol.add_col(COLM(edges[m]));
            thresh = std::max(thresh, edges[m].deathTime());
            continue;
        }
        tempCol.clear();
        return INFINITY;
    }
    tempCol.clear();
    return thresh;
}

// Squared diameter of the minimal enclosing ball of a triangle with squared side lengths x, y, z.
inline double triGrade(double x, double y, double z){
    const double m = std::max(x, std::max(y, z));
    if(2*m >= x + y + z) return m;   // right or obtuse: the longest side is a diameter
    return 2*x*y*z/(x*y + y*z + z*x - (x*x + y*y + z*z)/2);
}

// Skipping cofaces: since the boundary of the boundary of abcd is zero, the boundary of abd is the sum of
// the boundaries of abc, acd and bcd. If those three triangles all have strictly lower grade than abd,
// then removing abd from the filtration does not change the space of boundaries at any grade, so the H1
// barcode is unchanged, and abd does not need to be reduced. (Strictly lower grade makes the removals
// well-founded: a skipped triangle is never needed to justify skipping a triangle of lower or equal grade.)
// While going through the cofaces of an edge ab, the points c for which abc was reduced to zero are kept
// as candidates. To be safe against rounding (the grade of the same triangle computed from different
// edges can differ in the last bits), "strictly lower" is checked with a small relative margin.
const double skipMargin = 1 - 1e-12;
long long noOfSkipped = 0, noOfReduced = 0;

// While sorting, the squared length of an edge is kept in the first 8 bytes (edgeInfo1 and edgeInfo2, unused at that time).
inline void setSortKey(edge& e, double d){
    uint64_t b; std::memcpy(&b, &d, sizeof(double));
    e.edgeInfo1 = (int32_t)(uint32_t)b;
    e.edgeInfo2 = (int32_t)(uint32_t)(b >> 32);
}
inline double sortKey(const edge& e){
    const uint64_t b = (uint64_t)(uint32_t)e.edgeInfo1 | ((uint64_t)(uint32_t)e.edgeInfo2 << 32);
    double d; std::memcpy(&d, &b, sizeof(double)); return d;
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
            const double d = sqDist(i, j);
            setSortKey(e, d);
            e.slot = -1;
            e.points = {(uint16_t)i, (uint16_t)j};
            ijCounter ++;
        }
    }
    auto byDist = [](const edge& lhs, const edge& rhs) { return sortKey(lhs) < sortKey(rhs); };
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
                    edg.slot = (int32_t)deathTimes.size();
                    deathTimes.push_back(-1);
                    slotEdge.push_back((int32_t)k);

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
    struct Cand { int32_t c; double g, pc, qc; };   // c, grade of pqc, |pc|^2, |qc|^2
    std::vector<Cand> zeroCofaces;   // candidates c: cofaces pqc of the current edge pq that reduced to zero
    index i = -1;
    for(edge& edg : edges){
        i++;
        double pq = edg.distance();
        bitColumnSingleEdge.add_index(i);
        if(reduceUnc(&bitColumnSingleEdge, pq, true, 1)){
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
        zeroCofaces.clear();
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
            // Skip pqr if some earlier zero coface pqc has pqc, prc and qrc all of strictly lower grade (see skipMargin).
            // The most recent candidates are tried first, since they are the most likely to work.
            {
                const double gMax = pqrDiameter*skipMargin;
                bool skip = false;
                for(int t = (int)zeroCofaces.size() - 1; t >= 0; t--){
                    const Cand& cand = zeroCofaces[t];
                    if(!(cand.g < gMax)) continue;
                    const double cr = sqDist(cand.c, rIndex);
                    if(triGrade(cand.pc, cr, pr) < gMax && triGrade(cand.qc, cr, qr) < gMax){
                        skip = true;
                        break;
                    }
                }
                if(skip){
                    noOfSkipped++;
                    continue;
                }
                noOfReduced++;
            }
            index prIndex = edgeIndices[pairIndex(pIndex, rIndex)];
            index qrIndex = edgeIndices[pairIndex(qIndex, rIndex)];
            bitColumnpqr.add_index(prIndex);
            bitColumnpqr.add_index(qrIndex);
            bitColumnpqr.add_col(columnSingleEdge);
            const bool changedSomething = reduceUnc(&bitColumnpqr, pqrDiameter, false, 0);
            if(!changedSomething){
                zeroCofaces.push_back({(int32_t)rIndex, pqrDiameter, pr, qr});
            }
            if(changedSomething){
                threshold = threshForZeroingOut(i);
                if(threshold < diamThreshold - 0.00001){
                    break;
                }
            }
        }
    }
    for(index i=0; i < H1pivotcounter; i++){
        barcode[i] = {edges[slotEdge[i]].distance(), deathTimes[i]};
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
    exhCol.init(noOfPairs);
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
    std::cerr << "Cofaces skipped: " << noOfSkipped << ", reduced: " << noOfReduced << "\n";

    // Bars are (birth^2, death^2); take square roots for the actual barcode.
    // A death of -1 (never killed) is kept negative, so write_barcode still writes it as 'inf'.
    for(auto& [birth, death] : barcode){
        birth = std::sqrt(birth);
        if(death >= 0) death = std::sqrt(death);
    }
    if(!write_barcode(output_filename, barcode)){
        return 1;
    }
    return 0;
}