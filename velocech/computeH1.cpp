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

std::vector<std::vector<double>> pointCloud;
std::vector<std::pair<double, std::pair<index, index>>> sortedDistances;    //
std::vector<std::vector<std::pair<double, index>>> distancesFromEachP;
std::vector<std::pair<index, index>> edgeInfo;  //
std::vector<std::vector<index>> edgeIndices;
std::vector<column> columns;                    //
std::vector<double> deathTime;                  //
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
column testCol;

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

// Function to calculate the distance between two points. Can probably be deleted
double distpq(const std::vector<double>* p1, const std::vector<double>* p2) {
    int dim = (*p1).size();
    double distanceSquared = 0;
    for (int l=0; l<dim; l++){
        double diff = (*p1)[l] - (*p2)[l];
        distanceSquared += diff*diff;
    }
    return distanceSquared;
}
/*
double pqrDiameter(const std::vector<double>* p, const std::vector<double>* q, const std::vector<double>* r, double pq, int dim){
    double pr = distpq(p, r);
    double qr = distpq(q, r);
    return 2*pq*pr*qr/(-(pq*pq+pr*pr+qr*qr)/2 + pq*pr+pr*qr+qr*pq);
}*/

// can probably delete
void addprq(index pIndex, index rIndex, index qIndex, bit_tree_column* bitcol){
    (*bitcol).add_index(edgeIndices[pIndex][rIndex]);
    (*bitcol).add_index(edgeIndices[qIndex][rIndex]);
}

std::pair<index,index> prqPair(index pIndex, index rIndex, index qIndex){
    index prIndex = edgeIndices[pIndex][rIndex];
    index qrIndex = edgeIndices[qIndex][rIndex];
    return {std::min(prIndex, qrIndex), std::max(prIndex, qrIndex)};
}

// input the pointCloud indices instead?
bool obtuse(const std::vector<double>* p, const std::vector<double>* q, const std::vector<double>* r){
    double dotProduct1 = 0;
    double dotProduct2 = 0;
    double dotProduct3 = 0;
    //use distances like in computemu instead?
    for (int d=0; d < dim; d++){
        dotProduct1 += ((*q)[d] - (*p)[d])*((*r)[d] - (*p)[d]);
        dotProduct2 += ((*p)[d] - (*q)[d])*((*r)[d] - (*q)[d]);
        dotProduct3 += ((*q)[d] - (*r)[d])*((*p)[d] - (*r)[d]);
    }
    if(dotProduct1 <= 0 || dotProduct2 <= 0 || dotProduct3 <= 0){
        return true;
    }
    return false;
}

//returns true iff edge i is important
bool computemu(const index i){
    double distpq = sortedDistances[i].first;
    auto [pIndex, qIndex] = sortedDistances[i].second;
    for (index j=0; j < noOfPoints-1; j++){
        index m = distancesFromEachP[pIndex][j].second;
        if (m == qIndex) {
            //There is no mu, so ab is important
            edgeInfo[i].first = -3;
            return true;
        }
        /*
        double dotProduct = 0;
        for (int d=0; d < dim; d++){
            dotProduct += (pointCloud[pIndex][d] - pointCloud[m][d])*(pointCloud[qIndex][d] - pointCloud[m][d]);
        }*/
        index mpIndex = edgeIndices[m][pIndex];
        index mqIndex = edgeIndices[m][qIndex];
        //think about what to do in non-generic cases, e.g., four points forming a square
        // 
        if (sortedDistances[mpIndex].first + sortedDistances[mqIndex].first < distpq) {
            //not important. Not ordered atm.
            edgeInfo[i].first = mpIndex;
            edgeInfo[i].second = mqIndex;
            //edgeInfo[i] = {mpIndex, mqIndex};
            return false;
        }
    }
    return true;
}

// should only be called if a >= -1, which is the same as i not being important
bool simplifyNonImp(index i){
    auto [a, b] = edgeInfo[i];
    bool madeAChange = false;
    if(a >= 0){
        if(edgeInfo[a].first == -1){
            edgeInfo[i].first = edgeInfo[a].second;
            madeAChange = true;
        }
        else if(edgeInfo[a].first == -2){
            if(deathTime[a] <= std::max(sortedDistances[i].first, lowerThresh)){
                if(columns[a].size() == 0){
                    //std::cout << "kejggejo" << std::endl;
                    edgeInfo[i].first = -1;
                    madeAChange = true;
                }
                if(columns[a].size() == 1){
                    edgeInfo[i].first = columns[a][0];
                    madeAChange = true;
                }
            }
        }
        else if(edgeInfo[a].second == -1){
            edgeInfo[i].first = edgeInfo[a].first;
            madeAChange = true;
        }
    }
    if(b >= 0){
        if(edgeInfo[b].first == -1){
            edgeInfo[i].second = edgeInfo[b].second;
            madeAChange = true;
        }
        else if(edgeInfo[b].first == -2){
            if(deathTime[b] <= std::max(sortedDistances[i].first, lowerThresh)){
                if(columns[b].size() == 0){
                    edgeInfo[i].second = -1;
                    madeAChange = true;
                }
                if(columns[b].size() == 1){
                    edgeInfo[i].second = columns[b][0];
                    madeAChange = true;
                }
            }
        }
        else if(edgeInfo[b].second == -1){
            edgeInfo[i].second = edgeInfo[b].first;
            madeAChange = true;
        }
    }
    /*
    auto [x, y] = edgeInfo[i];
    if(x == -1 && y >= 0 && edgeInfo[y].first >= -1){
        edgeInfo[i] = edgeInfo[y];
        madeAChange = true;
    }
    else if(y == -1 && x >= 0 && edgeInfo[x].first >= -1){
        edgeInfo[i] = edgeInfo[x];
        madeAChange = true;
    }*/
    return madeAChange;
}

// should only be called for a == -2, which means i is important and its column has been initialized
// assumes columns are stored without the max el.
void simplifyImp(index i){
    if(columns[i].empty()){
        return;
    }
    index m = columns[i].back();
    // b is dead, or is important but not initialized, or
    // column of b is to the right of that of i, so we cannot add it.
    // make sure that deathTime[i] has gotten a value
    if(edgeInfo[m].first == -3 || deathTime[m] > std::max(deathTime[i], lowerThresh)){
        //exit
        return;
    }
    columns[i].pop_back();
    if(edgeInfo[m].first == -1 && edgeInfo[m].second == -1){
        //exit
        return;
    }
    //(*bitCol).clear();       // delete this if I know bitCol is empty at the beginning
    tempCol.add_col(columns[i]);
    if(edgeInfo[m].first == -2){
        //(*bitCol).add_col(columns[(*bitCol).get_max_index()]);
        tempCol.add_col(columns[m]);
    }
    else{
        if(edgeInfo[m].first >= 0){
            tempCol.add_index(edgeInfo[m].first);
        }
        if(edgeInfo[m].second >= 0){
            tempCol.add_index(edgeInfo[m].second);
        }
    }
    tempCol.get_col_and_clear(columns[i]);
}

// should only be called for a >= -1
void reduceNonImp(index i){
    //auto [a, b] = edgeInfo[i];
    while(simplifyNonImp(i)){}
}

// do I remove the max el here?
// returns true if the bitcol is 0 in the end
bool reduceImp(bit_tree_column* bitCol, double& deathT, bool doNotTouchOtherColumns, int inst){
    while(!bitCol->is_empty()){
        switch(inst){
            case 0:
                globalCounter++;
                break;
            case 1:
                globalCounter1++;
                break;
            case 2:
                globalCounter2++;
                break;
            case 3:
                globalCounter3++;
        }
        index m = bitCol->get_max_index();
        // m not important
        if(edgeInfo[m].first >= -1){
            bitCol->add_index(m);
            simplifyNonImp(m);
            if(edgeInfo[m].first >= 0){
                bitCol->add_index(edgeInfo[m].first);
            }
            if(edgeInfo[m].second >= 0){
                bitCol->add_index(edgeInfo[m].second);
            }
            continue;
        }
        // add column of m normally
        if(edgeInfo[m].first == -2 && deathTime[m] <= deathT){
            //globalCounter2 ++;
            bitCol->add_index(m);
            simplifyImp(m);
            bitCol->add_col(columns[m]);
            continue;
        }
        if(doNotTouchOtherColumns){
            return false;
        }
        bitCol->add_index(m);
        // found new pivot
        if(edgeInfo[m].first == -3){
            deathTime[m] = deathT;
            (*bitCol).get_col_and_clear(columns[m]);
            edgeInfo[m].first = -2;
            //std::cout << columns[m].size() << std::endl;
            return true;
        }
        // swap with column that has later deathTime
        // initialize deathTimes as -1?
        // maybe this can be sped up by not copying a whole bitcolumn
        if(deathTime[m] > deathT){
            std::swap(deathTime[m], deathT);
            // const double temp = deathTime[m];
            // deathTime[m] = deathT;
            // deathT = temp;
            column tempcol;
            bitCol->get_col(tempcol);
            bitCol->add_col(columns[m]);
            std::swap(columns[m], tempcol);
            //tempCol = (*bitCol);
            //bitCol->add_col(columns[m]);
            //tempCol.get_col_and_clear(columns[m]);
        }
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
        if(edgeInfo[m].first >= -1){
            tempCol.add_index(m);
            if(edgeInfo[m].first >= 0){
                tempCol.add_index(edgeInfo[m].first);
            }
            if(edgeInfo[m].second >= 0){
                tempCol.add_index(edgeInfo[m].second);
            }
            continue;
        }
        // add column of m normally
        if(edgeInfo[m].first == -2){
            //globalCounter2 ++;
            tempCol.add_index(m);
            tempCol.add_col(columns[m]);
            thresh = std::max(thresh, deathTime[m]);
            continue;
        }
        tempCol.clear();
        return INFINITY;
    }
    tempCol.clear();
    return thresh;
}

//One long vector with an element for each pair i<j of point indices, sorted by sq distance.  , std::vector<index>* edgeIndices
void allSortedDistances() {
    auto startx = std::chrono::steady_clock::now();
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
    for (index k=0; k < noOfPairs; k++){
        index i = sortedDistances[k].second.first;
        index j = sortedDistances[k].second.second;
        edgeIndices[i][j] = k;
        edgeIndices[j][i] = k;
    }
    //std::cout << sortedDistances[0].second.first << " here2" << std::endl;
    auto stopx = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_secondsx = stopx - startx;
    std::cout << "allSortedDistances: " << elapsed_secondsx.count() << " seconds\n";
}

// sortedDistancesFromEachPoint2 in test.cpp
std::vector<std::vector<std::pair<double, index>>> sortedDistancesFromEachPoint() {
    auto startx = std::chrono::steady_clock::now();
    std::vector<std::vector<std::pair<double, index>>> distFromEachPoint(noOfPoints);
    for (index i=0; i < noOfPairs; i++){
        double dist = sortedDistances[i].first;
        auto [a, b] = sortedDistances[i].second;
        distFromEachPoint[a].push_back({dist,b});
        distFromEachPoint[b].push_back({dist,a});
    }
    auto stopx = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_secondsx = stopx - startx;
    std::cout << "sortedDistancesFromEachPoint: " << elapsed_secondsx.count() << " seconds\n";
    return distFromEachPoint;
}

std::vector<std::pair<double, double>> H1barcode(){
    std::vector<bool> dead(noOfPairs, false);
    std::vector<index> mu(noOfPairs, -1);
    std::vector<std::pair<index, index>> columnsH0(noOfPoints, std::pair{-1, -1});
    //std::vector<std::pair<index, index>> output(noOfPoints, std::pair{-1, -1});
    std::vector<index> H1PivotEdges(0);
    //consider not having allSortedDistances output distances, since I don't use them here
    auto H0startx = std::chrono::steady_clock::now();
    index jCounter = 0;
    // iterate over the elements instead of index?
    for(index i=0; i < noOfPairs; i++){
        if(!computemu(i)){
            // edge i is not important
            continue;
        }
        int pivot = sortedDistances[i].second.second;
        int nonpivot = sortedDistances[i].second.first;
        while(true){
            //If the pivot is already a pivot for an earlier column
            if(columnsH0[pivot].second != -1){
                int otherNonpivot = columnsH0[pivot].first;
                if(nonpivot == otherNonpivot){
                    //the column goes to zero; this means edge #i contributes to H1
                    H1pivotcounter++;
                    H1PivotEdges.push_back(i);
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
                edgeInfo[i].first = -1;
                edgeInfo[i].second = -1;
                //pivotsFound ++;
                break;
            }
        }
    }
    auto H0stopx = std::chrono::steady_clock::now();
    std::chrono::duration<double> H0elapsed_secondsx = H0stopx - H0startx;
    std::cout << "computing H0: " << H0elapsed_secondsx.count() << " seconds\n";
    std::cout << "jCounter: " << jCounter << std::endl;
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
    for(index i=0; i < noOfPairs; i++){
        if(edgeInfo[i].first >= -1){
            //reduceNonImp(i);
            while(simplifyNonImp(i)){}
            //std::cout << edgeInfo[i].first << ", " << edgeInfo[i].second << std::endl;
            if((edgeInfo[i].first == -1) && (edgeInfo[i].second == -1)){
                zeroCounter ++;
                continue;
            }
        }
        double pq = sortedDistances[i].first;
        lowerThresh = pq;
        bitColumnSingleEdge.add_index(i);
        if(reduceImp(&bitColumnSingleEdge, pq, true, 1)){
            //zeroCounter ++;
            edgeInfo[i].first = -1;
            edgeInfo[i].second = -1;
            continue;
        }
        const index pIndex = sortedDistances[i].second.first;
        const index qIndex = sortedDistances[i].second.second;
        const index mIndex = mu[i];
        const std::vector<double> p = pointCloud[pIndex];
        const std::vector<double> q = pointCloud[qIndex];
        column columnSingleEdge = {};
        //index mxx = bitColumnSingleEdge.get_max_index();
        //threshRedCol = bitColumnSingleEdge;
        index redmax = bitColumnSingleEdge.get_max_index();
        bool importnt = edgeInfo[i].first <= -2;
        bitColumnSingleEdge.get_col_and_clear(columnSingleEdge);
        //go through the cofaces of edge i
        //------------------------------------
        //------------------------------------
        secondLoop ++;
        //column acute(0);
        double threshold = INFINITY;
        double diamThreshold = 0;
        int threshCounter = 0;
        //int threshCounter2 = 5;
        int counterX = 0;
        for(const std::pair<double, index>& rPair : distancesFromEachP[pIndex]){
            counterX ++;
            auto [pr, rIndex] = rPair;
            //if(threshCounter % 10 == 9 || counterX % 20 == 0){
            if(threshCounter % 10 == 9){
                double zeroThresh = threshForZeroingOut(i);
                //std::cout << threshCounter << ", " << zeroThresh << std::endl;
                if(zeroThresh < diamThreshold - 0.00001){
                    //std::cout << i << " and " << diamThreshold << " and " << zeroThresh << std::endl;
                    break;
                }
                //if(threshCounter % 10 == 9){threshCounter ++;}
                threshCounter ++;
                //threshCounter2 = threshCounter2*2;
            }
            if(threshCounter > 20 && pr + pq < diamThreshold - 0.00001){
                oblCounter1 ++;
                continue;
            }
            if(rIndex == qIndex){
                continue;
            }
            auto [prIndex, qrIndex] = prqPair(pIndex, rIndex, qIndex);
            const double qr = sortedDistances[qrIndex].first;
            // If the triangle is obtuse, it's not important, so we can skip it
            if(pq + pr <= qr || pq + qr <= pr){
                //oblCounter ++;
                continue;
            }
            // this is the square of the diameter of the circumcircle of pqr
            double pqrDiameter = 2*pq*pr*qr/(-(pq*pq+pr*pr+qr*qr)/2 + pq*pr+pr*qr+qr*pq);
            if(pr + qr <= pq){
                diamThreshold = std::max(diamThreshold, pqrDiameter);
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
            if(importnt){
                counterImportant ++;
            }
            if(pr > threshold){
                //tripleCounter++;
                break;
            }
            /*if(reduceImp(&threshRedCol, pr, true, 2)){
                counterImportant ++;
                break;
            }*/
            if(pqrDiameter > threshold){
                continue;
            }
            counter ++;
            // first check if only pq contributes to the column
            // bitColumnpqr is empty
            bitColumnpqr.add_index(prIndex);
            bitColumnpqr.add_index(qrIndex);
            // Only pq contributes, which means the cofaces of pq w/ diam>pqrDiameter will
            // have pq zero out
            if(reduceImp(&bitColumnpqr, pqrDiameter, true, 3)){
                threshold = pqrDiameter;
                additionCounter1 ++;
            }
            bitColumnpqr.add_col(columnSingleEdge);
            if(!reduceImp(&bitColumnpqr, pqrDiameter, false, 0)){
                std::cout << "bitColumnpqr is not empty, and it should be" << std::endl;
            }
        }
    }
    //std::cout << "Barcode with squared values:" << std::endl;
    for(index i=0; i < H1pivotcounter; i++){
        barcode[i] = {sortedDistances[H1PivotEdges[i]].first, deathTime[H1PivotEdges[i]]};
        //std::cout << barcode[i].first << ", " << barcode[i].second << std::endl;
        // if(barcode[i].second > 500){
        //     std::cout << barcode[i].first << ", " << barcode[i].second << std::endl;
        // }
        if(barcode[i].first > barcode[i].second){
            std::cout << "-------death before birth, oh no!---------" << std::endl;
        }
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
    //std::cout << "globalCounter2:   " << globalCounter2 << std::endl;
    std::cout << "globalCounter3:   " << globalCounter3 << std::endl;
    std::cout << "H1pivotcounter:   " << H1pivotcounter << std::endl;
    return barcode;
}


int main() {
    double min_coord = 0.0;
    double max_coord = 60.0;
    dim = 3;
    noOfPoints = 3000;

    std::vector<double> point;
    std::vector<double> center = {30, 30};

    noOfPairs = noOfPoints*(noOfPoints-1)/2;
    edgeIndices.resize(noOfPoints, std::vector<index>(noOfPoints));
    edgeInfo.resize(noOfPairs, {-3, -3});
    columns.resize(noOfPairs);
    deathTime.resize(noOfPairs, -1);
    tempCol.init(noOfPairs);

    pointCloud.resize(noOfPoints, std::vector<double>(dim));
    for (int i=0; i < noOfPoints; i++){
        while(true){
            point = generateRandomPoint(min_coord, max_coord);
            std::vector<double> twoD_proj_of_point = point;
            twoD_proj_of_point.resize(2);
            if(440 > distpq(&twoD_proj_of_point, &center) && distpq(&twoD_proj_of_point, &center) > 400){
            //if(true){
                pointCloud[i] = point;
                //std::cout << point[0] << ", " << point[1] << ", " << distpq(&point, &center) << std::endl;
                break;
            }
        }
    }

    allSortedDistances();
    distancesFromEachP = sortedDistancesFromEachPoint();
    //std::vector<index> edgeIndices(noOfPairs);
    //std::vector<std::pair<double, std::pair<index, index>>> sD = allSortedDistances(&pointCloud, dimensions, &edgeIndices);
    //std::vector<std::pair<double, double>> barcode = H0barcode(&pointCloud, ...);
    //std::vector<std::pair<double, double>> barcode = 
    auto start = std::chrono::steady_clock::now();
    H1barcode();
    //std::vector<std::pair<int, int>> H0 = H0pivots(&pointCloud);
    //std::vector<std::vector<std::pair<double, int>>> asd = allSortedDistances(&pointCloud, dimensions);
    //std::vector<std::pair<double, std::pair<int, int>>> asd2 = allSortedDistancesGlobal(&pointCloud, dimensions);
    //std::vector<std::vector<std::pair<double, int>>> asd = firstNSortedDistances(&pointCloud, 10);
    /*
    std::cout << "---------"<< endl;
    for (int i=0; i < barcode.size(); i++){
        //for (int j=0; j < noOfPoints-1; j++){
        std::cout << barcode[i].first << ", " << barcode[i].second << std::endl;
        //}
    }
    */
    //auto mid1 = std::chrono::steady_clock::now();
    //vector<vector<int>> m = mu(&pointCloud, &asd);
    //auto mid2 = std::chrono::steady_clock::now();
    //vector<vector<vector<pair<int, int>>>> phi = phiMatrix(&m);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> totalTime = end - start;
    std::cout << "Total time: " << totalTime.count() << " seconds\n";
    return 0;
}