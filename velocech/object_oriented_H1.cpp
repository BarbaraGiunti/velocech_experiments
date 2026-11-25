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
using namespace phat;

class edge;

std::vector<std::vector<double>> pointCloud;
std::vector<std::pair<double, std::pair<index, index>>> sortedDistances;
std::vector<std::vector<std::pair<double, index>>> distancesFromEachP;
std::vector<std::vector<index>> edgeIndices;
int dim;
int32_t noOfPoints;
index noOfPairs;
bit_tree_column tempCol;
index H1pivotcounter = 0;
//
index globalCounter = 0;
index globalCounter1 = 0;
index globalCounter2 = 0;
index globalCounter3 = 0;
double lowerThresh = 0;

std::vector<edge> edges;

class edge
{
    protected:
    public:
        std::pair<index, index> points;
        double distance;
        column colm;
        double deathTime;
        index edgeInfo1;
        index edgeInfo2;
        void init(double dist, std::pair<index, index> ab){
            distance = dist;
            points = ab;
            deathTime = -1;
            edgeInfo1 = edgeInfo2 = -3;
        }
        
        // should only be called if edgeInfo1 >= -1, which is the same as i not being important
        void simplifyNonImp(){
            bool madeAChange = true;
            while(madeAChange){
                madeAChange = false;
                edge& edgeA = edges[edgeInfo1];
                if(edgeInfo1 >= 0 && edgeA.deathTime <= std::max(distance, lowerThresh)){
                    if(edgeA.edgeInfo1 == -1){
                        edgeInfo1 = edgeA.edgeInfo2;
                        madeAChange = true;
                    }
                    else if(edgeA.edgeInfo1 == -2){
                        if(edgeA.colm.size() == 0){
                            edgeInfo1 = -1;
                            madeAChange = true;
                        }
                        if(edgeA.colm.size() == 1){
                            edgeInfo1 = edgeA.colm[0];
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
                edge& edgeB = edges[edgeInfo2];
                if(edgeInfo2 >= 0 && edgeB.deathTime <= std::max(distance, lowerThresh)){
                    if(edgeB.edgeInfo1 == -1){
                        edgeInfo2 = edgeB.edgeInfo2;
                        madeAChange = true;
                    }
                    else if(edgeB.edgeInfo1 == -2){
                        if(edgeB.colm.size() == 0){
                            edgeInfo2 = -1;
                            madeAChange = true;
                        }
                        if(edgeB.colm.size() == 1){
                            edgeInfo2 = edgeB.colm[0];
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
        bool computemu(){
            double distpq = distance;
            auto [pIndex, qIndex] = points;
            for (std::pair<double, phat::index>& mPair : distancesFromEachP[pIndex]){
                index m = mPair.second;
                if (m == qIndex) {
                    //There is no mu, so ab is important
                    edgeInfo1 = -3;
                    return true;
                }
                index mpIndex = edgeIndices[m][pIndex];
                index mqIndex = edgeIndices[m][qIndex];
                //think about what to do in non-generic cases, e.g., four points forming a square
                // 
                if (edges[mpIndex].distance + edges[mqIndex].distance < distpq) {
                    //not important. Not ordered atm.
                    edgeInfo1 = mpIndex;
                    edgeInfo2 = mqIndex;
                    //edgeInfo[i] = {mpIndex, mqIndex};
                    return false;
                }
            }
            return true;
        }

        // should only be called for a == -2, which means i is important and its column has been initialized
        // assumes columns are stored without the max el.
        void simplifyImp(){
            if(colm.empty()){
                return;
            }
            const edge& edgeM = edges[colm.back()];
            // b is dead, or is important but not initialized, or
            // column of b is to the right of that of i, so we cannot add it.
            // make sure that deathTime[i] has gotten a value
            if(edgeM.edgeInfo1 == -3 || edgeM.deathTime > std::max(deathTime, lowerThresh)){
                //exit
                return;
            }
            colm.pop_back();
            if(edgeM.edgeInfo1 == -1 && edgeM.edgeInfo2 == -1){
                //exit
                return;
            }
            //(*bitCol).clear();       // delete this if I know bitCol is empty at the beginning
            tempCol.add_col(colm);
            if(edgeM.edgeInfo1 == -2){
                //(*bitCol).add_col(columns[(*bitCol).get_max_index()]);
                tempCol.add_col(edgeM.colm);
            }
            else{
                if(edgeM.edgeInfo1 >= 0){
                    tempCol.add_index(edgeM.edgeInfo1);
                }
                if(edgeM.edgeInfo2 >= 0){
                    tempCol.add_index(edgeM.edgeInfo2);
                }
            }
            tempCol.get_col_and_clear(colm);
        }
};

std::vector<double> generateRandomPoint(double min_val, double max_val) {
    // Use a random device to seed the Mersenne Twister engine
    static std::random_device rd;
    static std::mt19937 gen(rd());
    // Define the uniform real distribution for the desired range
    std::uniform_real_distribution<double> distrib(min_val, max_val);

    std::vector<double> point;
    for (int i = 0; i < dim; ++i) {
        // Generate a random coordinate and add it to the point
        point.push_back(distrib(gen));
    }
    return point;
}

// Function to calculate the distance between two points.
double distpq(const std::vector<double>* p1, const std::vector<double>* p2) {
    int dim = (*p1).size();
    double distanceSquared = 0;
    for (int l=0; l<dim; l++){
        double diff = (*p1)[l] - (*p2)[l];
        distanceSquared += diff*diff;
    }
    return distanceSquared;
}

int findIndex(const std::vector<std::pair<double, index>>& myVec, double target)
{
    auto it = std::upper_bound(
        myVec.begin(),
        myVec.end(),
        target,
        [](double value, const std::pair<double, index>& p) {
            return value < p.first;
        }
    );

    if (it == myVec.begin()) {
        return -1;  // No element <= query
    }

    --it;
    return static_cast<int>(it - myVec.begin());
}

// do I remove the max el here?
// returns true if the bitcol is 0 in the end
// If doNotTouchOtherColumns, then return true if the bitcol is 0 in the end
// If !doNotTouchOtherColumns, then return true if some column has gotten a new deathTime
// (or whatever the precise statement is)
bool reduceImp(bit_tree_column* bitCol, double& deathT, bool doNotTouchOtherColumns, int inst){
    bool returnVal = false;
    while(!bitCol->is_empty()){
        switch(inst){
            case 0:
                globalCounter++;
                break;
            case 1:
                globalCounter1++;
                break;
            case 2:
                // globalCounter2++;
                break;
            case 3:
                globalCounter3++;
        }
        index m = bitCol->get_max_index();
        // m not important
        if(edges[m].edgeInfo1 >= -1 && edges[m].deathTime <= deathT){
            globalCounter2 ++;
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
            bitCol->add_col(edges[m].colm);
            continue;
        }
        if(doNotTouchOtherColumns){
            return false;
        }
        bitCol->add_index(m);
        // found new pivot
        if(edges[m].edgeInfo1 == -3){
            edges[m].deathTime = deathT;
            (*bitCol).get_col_and_clear(edges[m].colm);
            if(edges[m].colm.size() == 0){
                edges[m].edgeInfo1 = -1;
                edges[m].edgeInfo2 = -1;
            }
            else if(edges[m].colm.size() == 1){
                edges[m].edgeInfo1 = edges[m].colm[0];
                edges[m].edgeInfo2 = -1;
            }
            else if(edges[m].colm.size() == 2){
                edges[m].edgeInfo1 = edges[m].colm[0];
                edges[m].edgeInfo2 = edges[m].colm[1];
            }
            else{
                edges[m].edgeInfo1 = -2;
            }
            return true;

        }
        // swap with column that has later deathTime
        if(edges[m].deathTime > deathT){
            returnVal = true;
            //std::swap(deathTime[m], deathT);
            const double temp = edges[m].deathTime;
            edges[m].deathTime = deathT;
            deathT = temp;
            //
            column tempcol;
            bitCol->get_col(tempcol);
            if(edges[m].edgeInfo1 == -2){
                bitCol->add_col(edges[m].colm);
            }
            else{
                if(edges[m].edgeInfo1 >= 0){
                    bitCol->add_index(edges[m].edgeInfo1);
                }
                if(edges[m].edgeInfo2 >= 0){
                    bitCol->add_index(edges[m].edgeInfo2);
                }
            }
            if(tempcol.size() == 0){
                edges[m].edgeInfo1 = -1;
                edges[m].edgeInfo2 = -1;
            }
            else if(tempcol.size() == 1){
                edges[m].edgeInfo1 = tempcol[0];
                edges[m].edgeInfo2 = -1;
            }
            else if(tempcol.size() == 2){
                edges[m].edgeInfo1 = tempcol[0];
                edges[m].edgeInfo2 = tempcol[1];
            }
            else{
                edges[m].edgeInfo1 = -2;
                std::swap(edges[m].colm, tempcol);
            }
            //tempCol = (*bitCol);
            //bitCol->add_col(columns[m]);
            //tempCol.get_col_and_clear(columns[m]);
        }
    }
    if(!doNotTouchOtherColumns){
        return returnVal;
    }
    return true;
}

double threshForZeroingOut(index i){
    double thresh = 0;
    tempCol.clear();
    tempCol.add_index(i);
    while(!tempCol.is_empty()){
        index m = tempCol.get_max_index();
        // m not important
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
            //globalCounter2 ++;
            tempCol.add_index(m);
            tempCol.add_col(edges[m].colm);
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
    auto startx = std::chrono::steady_clock::now();
    //std::vector<std::pair<double, std::pair<index, index>>> sortedDistances(noOfPairs);
    edges.resize(noOfPairs);
    //std::vector<std::pair<double, std::pair<index, index>>> sortedDistances(noOfPairs);
    sortedDistances.resize(noOfPairs);
    index ijCounter=0;
    for (index i=0; i < noOfPoints; i++){
        for (index j=i+1; j < noOfPoints; j++){
            sortedDistances[ijCounter] = {distpq(&pointCloud[i], &pointCloud[j]), {i,j}};
            ijCounter ++;
        }
    }
    std::sort(sortedDistances.begin(), sortedDistances.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
    });
    auto startxxx = std::chrono::steady_clock::now();
    for (index k=0; k < noOfPairs; k++){
        const auto& [dist, ab] = sortedDistances[k];
        edges[k].init(dist, ab);
        edgeIndices[ab.first][ab.second] = k;
        edgeIndices[ab.second][ab.first] = k;
    }
    //std::cout << " here2" << std::endl;
    auto stopx = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_secondsx = stopx - startx;
    std::chrono::duration<double> elapsed_secondsxxx = stopx - startxxx;
    std::cout << "allSortedDistances: " << elapsed_secondsx.count() << " seconds\n";
    std::cout << "allSortedDistances, part with edges: " << elapsed_secondsxxx.count() << " seconds\n";
}

// sortedDistancesFromEachPoint2 in test.cpp
std::vector<std::vector<std::pair<double, index>>> sortedDistancesFromEachPoint() {
    auto startx = std::chrono::steady_clock::now();
    std::vector<std::vector<std::pair<double, index>>> distFromEachPoint(noOfPoints);
    for (const auto& [dist, ab] : sortedDistances){
        distFromEachPoint[ab.first].push_back({dist,ab.second});
        distFromEachPoint[ab.second].push_back({dist,ab.first});
    }
    auto stopx = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_secondsx = stopx - startx;
    std::cout << "sortedDistancesFromEachPoint: " << elapsed_secondsx.count() << " seconds\n";
    return distFromEachPoint;
}

std::vector<std::pair<double, double>> H1barcode(){
    std::vector<std::pair<index, index>> columnsH0(noOfPoints, std::pair{-1, -1});
    std::vector<edge*> H1PivotEdges(0);
    auto H0startx = std::chrono::steady_clock::now();
    for(edge& edg : edges){
        //edge edgeI = edges[i];
        if(!edg.computemu()){
            // edge i is not important
            continue;
        }
        int pivot = edg.points.second;
        int nonpivot = edg.points.first;
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
                //output[pivot] = sortedDistances[i].second;
                //the edge is dead (so can be ignored when computing H1)
                edg.edgeInfo1 = -1;
                edg.edgeInfo2 = -1;
                //pivotsFound ++;
                break;
            }
        }
    }
    auto H0stopx = std::chrono::steady_clock::now();
    std::chrono::duration<double> H0elapsed_secondsx = H0stopx - H0startx;
    std::cout << "computing H0: " << H0elapsed_secondsx.count() << " seconds\n";
    /////////////////////////////////////
    ////         Computing H1        ////
    /////////////////////////////////////
    std::cout << noOfPairs << std::endl;
    column zeroColumn = {};
    //std::vector<std::pair<column, double>> edgeInfo(noOfPairs, std::pair{zeroColumn, INFINITY});
    auto startx = std::chrono::steady_clock::now();
    //std::vector<index> H1PivotEdges(H1pivotcounter);
    index pivotCounter = 0;
    std::vector<std::pair<double, double>> barcode(H1pivotcounter);
    //return barcode;
    index counter = 0;
    index counterImportant = 0;
    index quickIterCounter = 0;
    index additionCounter1 = 0;
    index additionCounter2 = 0;
    index tripleCounter = 0;
    index zeroCounter = 0;
    index oblCounter = 0;
    index oblCounter1 = 0;
    index oblCounter2 = 0;
    index secondLoop = 0;
    index canSkipCounter = 0;
    bit_tree_column bitColumnSingleEdge, bitColumnpqr;
    //threshRedCol;
    bitColumnSingleEdge.init(noOfPairs);
    bitColumnpqr.init(noOfPairs);
    index i = -1;
    for(edge& edg : edges){
        i++;
        if(edg.edgeInfo1 >= -1){
            //reduceNonImp(i);
            edg.simplifyNonImp();
            //std::cout << edgeInfo[i].first << ", " << edgeInfo[i].second << std::endl;
            if((edg.edgeInfo1 == -1) && (edg.edgeInfo2 == -1) && edg.deathTime < 0){
                zeroCounter ++;
                continue;
            }
        }
        double pq = edg.distance;
        lowerThresh = pq;
        bitColumnSingleEdge.add_index(i);
        if(reduceImp(&bitColumnSingleEdge, pq, true, 1)){
            //zeroCounter ++;
            edg.edgeInfo1 = -1;
            edg.edgeInfo2 = -1;
            continue;
        }
        const index pIndex = edg.points.first;
        const index qIndex = edg.points.second;
        const std::vector<double> p = pointCloud[pIndex];
        const std::vector<double> q = pointCloud[qIndex];
        column columnSingleEdge = {};
        // bool importnt = edg.edgeInfo1 <= -2;
        bitColumnSingleEdge.get_col_and_clear(columnSingleEdge);
        secondLoop ++;
        //column acute(0);
        double threshold = threshForZeroingOut(i);
        double diamThreshold = 0;
        int threshCounter = 0;
        //int threshCounter2 = 5;
        int counterX = 0;
        ////////////////////////////////////////////
        ////  Go through the cofaces of edge i  ////
        ////////////////////////////////////////////
        for(index j=0; j < distancesFromEachP[pIndex].size(); j++){
        // for(const std::pair<double, index>& rPair : distancesFromEachP[pIndex]){
            counterX ++;
            tripleCounter ++;
            // auto [pr, rIndex] = rPair;
            auto [pr, rIndex] = distancesFromEachP[pIndex][j];
            if(rIndex == qIndex || (threshCounter > 15 && pr + pq < diamThreshold - 0.00001)){
                j = std::max(static_cast<int>(j), findIndex(distancesFromEachP[pIndex], diamThreshold-pq));
                continue;
            }
            //////////////////////////////////////////////////////////////////////////
            ////  A bunch of geometric checks that may allow skipping iterations  ////
            //////////////////////////////////////////////////////////////////////////
            index prIndex = edgeIndices[pIndex][rIndex];
            index qrIndex = edgeIndices[qIndex][rIndex];
            const double qr = edges[qrIndex].distance;
            // If the triangle is obtuse, it's not important, so we can skip it
            if(pq + pr <= qr || pq + qr <= pr){
                //oblCounter ++;
                continue;
            }
            // this is the square of the diameter of the circumcircle of pqr
            double pqrDiameter = 2*pq*pr*qr/(-(pq*pq+pr*pr+qr*qr)/2 + pq*pr+pr*qr+qr*pq);
            if(pr + qr <= pq){
                diamThreshold = std::max(diamThreshold, pqrDiameter);
                if(threshold < diamThreshold - 0.00001){
                    break;
                }
                oblCounter ++;
                threshCounter ++;
                continue;
            }
            if(pqrDiameter < diamThreshold - 0.00001){
                oblCounter2 ++;
                continue;
            }
            const std::vector<double> r = pointCloud[rIndex];
            canSkipCounter ++;
            // if(importnt){
            //     counterImportant ++;
            // }
            if(pr > threshold){
                //tripleCounter++;
                break;
            }
            if(pqrDiameter > threshold){
                continue;
            }
            ////////////////////////////////////////////////////////////
            ////  Done with geometry, the rest is column reduction  ////
            ////////////////////////////////////////////////////////////
            counter ++;
            // bitColumnpqr is empty
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
    //std::cout << "Barcode with squared values:" << std::endl;
    for(index i=0; i < H1pivotcounter; i++){
        barcode[i] = {(*H1PivotEdges[i]).distance, (*H1PivotEdges[i]).deathTime};
        //std::cout << barcode[i].first << ", " << barcode[i].second << std::endl;
        if(barcode[i].second > 500){
            std::cout << barcode[i].first << ", " << barcode[i].second << std::endl;
            // std::cout << (*H1PivotEdges[i]).edgeInfo1 << ", " << (*H1PivotEdges[i]).edgeInfo2 << std::endl;
        }
        // if(barcode[i].first > barcode[i].second){
        //     std::cout << "-------death before birth, oh no!---------" << std::endl;
        // }
    }
    auto stopx = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_secondsx = stopx - startx;
    std::cout << "H1 function: " << elapsed_secondsx.count() << " seconds\n";
    std::cout << "tripleCounter:    " << tripleCounter << std::endl;
    //std::cout << "quickIterCounter: " << quickIterCounter << std::endl;
    std::cout << "secondLoop:       " << secondLoop << std::endl;
    std::cout << "oblCounter:       " << oblCounter << std::endl;
    std::cout << "oblCounter1:      " << oblCounter1 << std::endl;
    std::cout << "oblCounter2:      " << oblCounter2 << std::endl;
    std::cout << "counter:          " << counter << std::endl;
    std::cout << "canSkipCounter:   " << canSkipCounter << std::endl;
    std::cout << "counterImportant: " << counterImportant << std::endl;
    std::cout << "zeroCounter:      " << zeroCounter << std::endl;
    std::cout << "additionCounter1: " << additionCounter1 << std::endl;
    std::cout << "globalCounter:    " << globalCounter << std::endl;
    std::cout << "globalCounter1:   " << globalCounter1 << std::endl;
    std::cout << "globalCounter2:   " << globalCounter2 << std::endl;
    std::cout << "globalCounter3:   " << globalCounter3 << std::endl;
    std::cout << "H1pivotcounter:   " << H1pivotcounter << std::endl;
    return barcode;
}


int main() {
    auto start = std::chrono::steady_clock::now();
    double min_coord = 0.0;
    double max_coord = 60.0;
    dim = 3;
    noOfPoints = 2000;

    std::vector<double> point;
    std::vector<double> center = {30, 30};

    noOfPairs = noOfPoints*(noOfPoints-1)/2;
    edgeIndices.resize(noOfPoints, std::vector<index>(noOfPoints));
    tempCol.init(noOfPairs);

    pointCloud.resize(noOfPoints, std::vector<double>(dim));
    for (int i=0; i < noOfPoints; i++){
        while(true){
            point = generateRandomPoint(min_coord, max_coord);
            std::vector<double> twoD_proj_of_point = point;
            twoD_proj_of_point.resize(2);
            if(500 > distpq(&twoD_proj_of_point, &center) && distpq(&twoD_proj_of_point, &center) > 400){
            // if(true){
                pointCloud[i] = point;
                //std::cout << point[0] << ", " << point[1] << ", " << distpq(&point, &center) << std::endl;
                break;
            }
        }
    }

    allSortedDistances();
    distancesFromEachP = sortedDistancesFromEachPoint();
    H1barcode();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> totalTime = end - start;
    std::cout << "Total time: " << totalTime.count() << " seconds\n";
    return 0;
}