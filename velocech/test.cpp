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
//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Min_sphere_of_spheres_d.h>
//#include <CGAL/Min_sphere_of_points_d_traits_3.h>
#include <array>
//using namespace std;
using namespace phat;

/*
typedef  CGAL::Simple_cartesian<double>                   K;
typedef  CGAL::Min_sphere_of_points_d_traits_3<K,double>  Traits;
typedef  CGAL::Min_sphere_of_spheres_d<Traits>            Min_sphere;
typedef  K::Point_3                                       Point;
*/

// Function to generate a vector of N random points in a specified range
std::vector<double> generateRandomPoint(double min_val, double max_val, int dimensions) {
    // Use a random device to seed the Mersenne Twister engine
    static std::random_device rd;
    static std::mt19937 gen(rd());
    // Define the uniform real distribution for the desired range
    std::uniform_real_distribution<double> distrib(min_val, max_val);

    std::vector<double> point;
    for (int i = 0; i < dimensions; ++i) {
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

double pqrDiameter(const std::vector<double>* p, const std::vector<double>* q, const std::vector<double>* r, double pq, int dim){
    double pr = distpq(p, r);
    double qr = distpq(q, r);
    return 2*pq*pr*qr/(-(pq*pq+pr*pr+qr*qr)/2 + pq*pr+pr*qr+qr*pq);
}

/*
bool inCircumsphere3D(const std::vector<double>* p, const std::vector<double>* q, const std::vector<double>* r, const std::vector<double>* x){
    std::array<Point, 3> P;
    P.at(0) = Point((*p)[0], (*p)[1], (*p)[2]);
    P.at(1) = Point((*q)[0], (*q)[1], (*q)[2]);
    P.at(2) = Point((*r)[0], (*r)[1], (*r)[2]);
    Min_sphere mc(P.begin(), P.end());
    Min_sphere::Cartesian_const_iterator ccib = mc.center_cartesian_begin(), ccie = mc.center_cartesian_end();
    double distSq = 0;
    int i = 0;
    for( ; ccib != ccie; ++ccib){
        distSq += (*ccib-(*x)[i])*(*ccib-(*x)[i]);
        i++;
    }
    return distSq <= mc.radius()*mc.radius();
} */


/*
bool inCircumsphere(const std::vector<double>* p, const std::vector<double>* q, const std::vector<double>* r, double pq, const std::vector<double>* x, int dim){
    const double radius = pqrDiameter(p, q, r, pq, dim)/4;
    std::vector<double> m(dim, 0);
    std::vector<double> v(dim, 0);
    std::vector<double> ptor(dim, 0);
    std::vector<double> ptoq(dim, 0);
    std::vector<double> mtox(dim, 0);
    double uDotpToq(0);
    //double sqrtpq = std::sqrt(pq);
    for (int i=0; i<dim; i++){
        m[i] = ((*p)[i] + (*q)[i])/2;
        ptor[i] = (*r)[i] - (*p)[i];
        ptoq[i] = (*q)[i] - (*p)[i];
        uDotpToq += ptoq[i]*ptor[i]/pq;
    }
    for (int i=0; i<dim; i++){
        mtox[i] = ptor[i] - uDotpToq*ptoq[i];
        x[i] = (m)[i] + ptor[i];
    }
    v = 
    const double center = ;
    return distpq(center, x) <= radius;
}*/

double lowerThresh(const std::vector<double>* p, const std::vector<double>* q, const std::vector<double>* r, double pq, int dim){
    double pr = distpq(p, r);
    double qr = distpq(q, r);
    return 2*pq*pr*qr/(-(pq*pq+pr*pr+qr*qr)/2 + pq*pr+pr*qr+qr*pq) - pq;
}

int obtuse(const std::vector<double>* p, const std::vector<double>* q, const std::vector<double>* r, int dim){
    double dotProduct1 = 0;
    double dotProduct2 = 0;
    double dotProduct3 = 0;
    for (int d=0; d < dim; d++){
        dotProduct1 += ((*q)[d] - (*p)[d])*((*r)[d] - (*p)[d]);
        dotProduct2 += ((*p)[d] - (*q)[d])*((*r)[d] - (*q)[d]);
        dotProduct3 += ((*q)[d] - (*r)[d])*((*p)[d] - (*r)[d]);
    }
    if(dotProduct1 <= 0){
        return 1;
    }
    if(dotProduct2 <= 0){
        return 2;
    }
    if(dotProduct3 <= 0){
        return 3;
    }
    return 0;
}

// Function to calculate the distance between two points. Can probably be deleted
double testFunction(const std::vector<double>* p1) {
    return (*p1)[0];
}

//Creates a vector of pairs (distance^2, point) to point #pointIndex for every other point,
//sorted by sq distances
//Haven't tested if this is correct
std::vector<std::pair<double, int>> sortedDistancesToAPoint(std::vector<std::vector<double>>* pointCloud, int pointIndex, int dim) {
    int numberOfPoints = static_cast<int>((*pointCloud).size());
    std::vector<std::pair<double, int>> distVector(numberOfPoints-1);
    int j=0;
    for (int i=0; i < numberOfPoints && i!= pointIndex; i++){
        double distanceSquared = 0;
        for (int l=0; l<dim; l++){
            double diff = (*pointCloud)[i][l] - (*pointCloud)[pointIndex][l];
            distanceSquared += diff*diff;
        }
        distVector[j] = {distanceSquared, i};
        j += 1;
    }
    std::sort(distVector.begin(), distVector.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
    });
    return distVector;
}

//Need to update with dimension as input
std::vector<std::vector<std::pair<double, int>>> firstNSortedDistances(std::vector<std::vector<double>>* pointCloud, int N) {
    int numberOfPoints = static_cast<int>((*pointCloud).size());
    std::vector<std::vector<std::pair<double, int>>> allDistances(numberOfPoints, std::vector<std::pair<double, int>>(numberOfPoints-1));
    int jj;
    double distance;
    for (int i=0; i < numberOfPoints; i++){
        int k = 0;
        double maxValue = std::numeric_limits<double>::infinity();
        for (int j=0; j < numberOfPoints; j++){
            double dx = (*pointCloud)[i][0] - (*pointCloud)[j][0];
            double dy = (*pointCloud)[i][1] - (*pointCloud)[j][1];
            double dz = (*pointCloud)[i][2] - (*pointCloud)[j][2];
            distance = std::hypot(dx, dy, dz);
            if (i < j) {
                jj = j-1;
            }
            else if (i > j) {
                jj = j;
            }
            else {
                continue;
            }
            if (k >= 2*N) {
                //the next lines assume N not too big
                std::sort(allDistances[i].begin(), allDistances[i].begin() + 2*N, [](const auto& lhs, const auto& rhs) {
                    return lhs.first < rhs.first;
                });
                maxValue = allDistances[i][N-1].first;
                k = N;
            }
            if (distance < maxValue) {
                allDistances[i][jj] = allDistances[i][k];
                allDistances[i][k] = {distance, j};
                k += 1;
            }
            else {
                allDistances[i][jj] = {distance, j};
            }
        }
        std::sort(allDistances[i].begin(), allDistances[i].begin() + k, [](const auto& lhs, const auto& rhs) {
            return lhs.first < rhs.first;
        });
        //allDistances[i].resize(N);
    }
    std::cout << "First " << N << " sorted distances, " << numberOfPoints << " points." << std::endl;
    return allDistances;
}

//Creates a vector of vectors. The outer vector has one vector for each point.
//Each inner vector has a pair (sq distance, point) for every other point, sorted by distances
//Haven't tested if this is correct
std::vector<std::vector<std::pair<double, index>>> sortedDistancesFromEachPoint(const std::vector<std::vector<double>>* pointCloud, const int dim) {
    auto startx = std::chrono::steady_clock::now();
    index numberOfPoints = static_cast<index>((*pointCloud).size());
    std::vector<std::vector<std::pair<double, index>>> allSortedDistances(numberOfPoints);
    std::vector<std::pair<double, index>> distVector(numberOfPoints-1);
    for (index k=0; k < numberOfPoints; k++){
        index j=0;
        for (index i=0; i < numberOfPoints; i++){
            if (i!= k) {
                double distanceSquared = 0;
                for (int l=0; l<dim; l++){
                    double diff = (*pointCloud)[i][l] - (*pointCloud)[k][l];
                    distanceSquared += diff*diff;
                }
                distVector[j] = {distanceSquared, i};
                j += 1;
            }
        }
        //Sorting is the current bottleneck
        std::sort(distVector.begin(), distVector.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.first < rhs.first;
        });
        allSortedDistances[k] = distVector;
    }
    std::cout << "All sorted distances, " << numberOfPoints << " points." << std::endl;
    auto stopx = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_secondsx = stopx - startx;
    std::cout << "sortedDistancesFromEachPoint: " << elapsed_secondsx.count() << " seconds\n";
    return allSortedDistances;
}

std::vector<std::vector<std::pair<double, index>>> sortedDistancesFromEachPoint2(const std::vector<std::vector<double>>* pointCloud, const int dim, std::vector<std::pair<double, std::pair<index, index>>>* sortedDistances) {
    auto startx = std::chrono::steady_clock::now();
    index numberOfPoints = static_cast<index>((*pointCloud).size());
    index noOfPairs = numberOfPoints*(numberOfPoints-1)/2;
    std::vector<std::vector<std::pair<double, index>>> distFromEachPoint(numberOfPoints);
    for (index i=0; i < noOfPairs; i++){
        double dist = (*sortedDistances)[i].first;
        index a = (*sortedDistances)[i].second.first;
        index b = (*sortedDistances)[i].second.second;
        // push_back is probably more expensive than initializing the vectors with the right
        // lengths. But this is not the bottleneck of the full algorithm.
        distFromEachPoint[a].push_back({dist,b});
        distFromEachPoint[b].push_back({dist,a});
    }
    auto stopx = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_secondsx = stopx - startx;
    std::cout << "sortedDistancesFromEachPoint2: " << elapsed_secondsx.count() << " seconds\n";
    return distFromEachPoint;
}

//One long vector with an element for each pair i<j of point indices, sorted by sq distance.  , std::vector<index>* edgeIndices
std::vector<std::pair<double, std::pair<index, index>>> allSortedDistances(const std::vector<std::vector<double>>* pointCloud, const int dim, std::vector<index>* edgeIndices) {
    auto startx = std::chrono::steady_clock::now();
    index numberOfPoints = static_cast<index>((*pointCloud).size());
    index numberOfUnorderedPairs = numberOfPoints*(numberOfPoints-1)/2;
    std::vector<std::pair<double, std::pair<index, index>>> sortedDistances(numberOfUnorderedPairs);
    // don't know how to enforce the right length on edgeIndices, so just assume for now it has the right length
    //std::vector<int>* edgeIndices(numberOfUnorderedPairs);
    //std::vector<std::pair<double, int>> distVector(numberOfPoints-1);
    index ijCounter=0;
    for (index i=0; i < numberOfPoints-1; i++){
        for (index j=i+1; j < numberOfPoints; j++){
            /* double distanceSquared = 0;
            for (int l=0; l<dim; l++){
                double diff = (*pointCloud)[i][l] - (*pointCloud)[j][l];
                distanceSquared += diff*diff;
            } */
            //double distanceSquared = distpq(&(*pointCloud)[i], &(*pointCloud)[j], dim);
            sortedDistances[ijCounter] = {distpq(&(*pointCloud)[i], &(*pointCloud)[j]), {i,j}};
            ijCounter ++;
        }
    }
    std::sort(sortedDistances.begin(), sortedDistances.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
    });
    for (index k=0; k < numberOfUnorderedPairs; k++){
        index i = sortedDistances[k].second.first;
        index j = sortedDistances[k].second.second;
        (*edgeIndices)[i + j*(j - 1)/2] = k;
    }
    //std::cout << sortedDistances[0].second.first << " here2" << std::endl;
    auto stopx = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_secondsx = stopx - startx;
    std::cout << "allSortedDistances: " << elapsed_secondsx.count() << " seconds\n";
    return sortedDistances;
}

/*
//not correct, I think, see comment inside. Might delete later idk.
vector<int> mu(vector<vector<double>>* pointCloud, vector<vector<pair<double, int>>>* distances, int dim) {
    int noOfPoints = static_cast<int>((*pointCloud).size());
    int noOfPairs = noOfPoints*(noOfPoints-1)/2;
    vector<int> mu(noOfPoints);
    vector<double> ijVector(dim);
    double dotProduct;
    int k = 0;
    for (int j=1; j < noOfPoints; j++){
        for (int i=0; i < j; i++){
            for (int d=0; d < dim; d++){
                ijVector[d] = (*pointCloud)[j][d] - (*pointCloud)[i][d];
            }
            for (int m=0; m < noOfPoints-1; m++){
                if ((*distances)[i][m].second == j) {
                    //There is no mu, so ij is important
                    //k can be computed as i + j(j-1)/2
                    mu[k] = -1;
                    k += 1;
                    break;
                }
                dotProduct = 0;
                //I think the following is wrong, m should be (*distances)[i][m].second or smth
                for (int d=0; d < dim; d++){
                    dotProduct += ijVector[d]*((*pointCloud)[m][d] - (*pointCloud)[i][d]);
                }
                if (dotProduct >= 0) {
                    //The dot product ij*im is nonnegative, so m is in the closed circumdisk of i, j.
                    //k can be computed as i + j(j-1)/2
                    mu[k] = (*distances)[i][m].second;
                    k += 1;
                    break;
                }
            }
        }
    }
    return mu;
}

vector<pair<int, int>> phi(vector<vector<int>>* mu, vector<vector<vector<pair<int, int>>>>* phiMatrix, int i, int j){
    if ((*phiMatrix)[i][j].size()) {
        return (*phiMatrix)[i][j];
    }
    int m = (*mu)[i][j];
    if (m == -1) {
        (*phiMatrix)[i][j] = (*phiMatrix)[j][i] = {{i,j}};
        return {{i,j}};
    }
    vector<pair<int, int>> phi_of_i_m = phi(mu, phiMatrix, i, m);
    vector<pair<int, int>> phi_of_m_j = phi(mu, phiMatrix, m, j);
    phi_of_i_m.insert(phi_of_i_m.end(), phi_of_m_j.begin(), phi_of_m_j.end());
    (*phiMatrix)[i][j] = (*phiMatrix)[j][i] = phi_of_i_m;
    return phi_of_i_m;
    return {{i,j}};
}

vector<vector<vector<pair<int, int>>>> phiMatrix(vector<vector<int>>* mu){
    int numberOfPoints = static_cast<int>((*mu).size());
    vector<vector<vector<pair<int, int>>>> phiMatrix(numberOfPoints, vector<vector<pair<int, int>>>(numberOfPoints, vector<pair<int, int>>(0)));
    for (int i=0; i<numberOfPoints; i++) {
        for (int j=i+1; j<numberOfPoints; j++) {
            if (!phiMatrix[i][j].size()) {
                phi(mu, &phiMatrix, i, j);
            }
        }
    }
    return phiMatrix;
}
*/
//std::vector<std::pair<int, int>>
int H0pivots(const std::vector<std::vector<double>>* pointCloud, std::vector<std::pair<column, double>>* edgeInfo, std::vector<std::vector<std::pair<double, index>>>* distancesFromEachP, std::vector<std::pair<double, std::pair<index, index>>>* sortedDistances, std::vector<index>* edgeIndices){
    const int dim = static_cast<int>((*pointCloud)[0].size());
    const index noOfPoints = static_cast<index>((*pointCloud).size());
    const index noOfPairs = noOfPoints*(noOfPoints-1)/2;
    //std::vector<std::pair<double, std::pair<index, index>>> sD = allSortedDistances(pointCloud, dim, edgeIndices);
    *sortedDistances = allSortedDistances(pointCloud, dim, edgeIndices);
    //std::cout << sD[0].second.first << std::endl << std::endl << std::endl;
    //*sortedDistances = sD;
    //std::cout << (*sortedDistances)[0].second.first << std::endl << std::endl << std::endl;
    std::vector<std::pair<index, index>> columns(noOfPoints, std::pair{-1, -1});
    std::vector<std::pair<index, index>> output(noOfPoints, std::pair{-1, -1});
    //vector with an entry for every edge, sorted by length.
    //bool = true if not dead. int refers to mu if negative, points to entry of H1edges if nonneg.
    //This is now given as input. How do I handle it here? atm declaring it in the H1 function
    //std::vector<std::pair<bool, int>> (*edgeInfo)(noOfPairs, std::pair{true, 0});
    //int pivotsFound = 0;
    //consider not having allSortedDistances output distances, since I don't use them here
    //*distancesFromEachP = sortedDistancesFromEachPoint(pointCloud, dim);
    *distancesFromEachP = sortedDistancesFromEachPoint2(pointCloud, dim, sortedDistances);
    auto startx = std::chrono::steady_clock::now();
    std::vector<index> testing(noOfPairs, -1);
    std::vector<std::pair<index, index>> testing2(noOfPairs);
    // currently outputting H1pivotcounter
    index H1pivotcounter = 0;
    index jCounter = 0;
    for(index i=0; i < noOfPairs; i++){
        const index pIndex = (*sortedDistances)[i].second.first;
        const index qIndex = (*sortedDistances)[i].second.second;
        bool important = false;
        for (index j=0; j < noOfPoints-1; j++){
            jCounter++;
            index m = (*distancesFromEachP)[pIndex][j].second;
            if (m == qIndex) {
                //There is no mu, so ab is important
                important = true;
                break;
            }
            double dotProduct = 0;
            for (int d=0; d < dim; d++){
                dotProduct += ((*pointCloud)[pIndex][d] - (*pointCloud)[m][d])*((*pointCloud)[qIndex][d] - (*pointCloud)[m][d]);
            }
            if (dotProduct <= 0) {
                testing2[i].first = i;
                testing2[i].second = 100;
                //The dot product mp*mq is nonpositive, so m is in the closed circumdisk of p, q.
                // we choose m = mu(p,q) and let (*edgeInfo)[i].first be pq+pm+qm
                // assuming that triangleColumn is set to zero and that add_index is ok with int input
                //bit_tree_column triangleColumn;
                //triangleColumn.add_index(i);
                index smallpm = std::min(pIndex, m);
                index largepm = std::max(pIndex, m);
                index smallqm = std::min(qIndex, m);
                index largeqm = std::max(qIndex, m);
                index indexpm = (*edgeIndices)[smallpm + largepm*(largepm - 1)/2];
                //triangleColumn.add_index((*edgeIndices)[indexpm]);
                index indexqm = (*edgeIndices)[smallqm + largeqm*(largeqm - 1)/2];
                break;
                //triangleColumn.add_index((*edgeIndices)[indexqm]);
                (*edgeInfo)[i].first = {std::min(indexpm, indexqm), std::max(indexpm, indexqm), i};
                (*edgeInfo)[i].second = (*sortedDistances)[i].first;
                /* double distpm = distpq((*pointCloud)[pIndex], (*pointCloud)[m], dim);
                double distqm = distpq((*pointCloud)[qIndex], (*pointCloud)[m], dim);
                double distpq = distpq((*pointCloud)[pIndex], (*pointCloud)[qIndex], dim);
                if(indexpm > i || indexqm > i){
                    std::cout << "--------------" << std::endl;
                    std::cout << "--------------" << std::endl;
                } */
                break;
            }
        }
        if(!important){
            continue;
        }
        int pivot = (*sortedDistances)[i].second.second;
        int nonpivot = (*sortedDistances)[i].second.first;
        while(true){
            //If the pivot is already a pivot for an earlier column
            if(columns[pivot].second != -1){
                int otherNonpivot = columns[pivot].first;
                if(nonpivot == otherNonpivot){
                    //the column goes to zero; this means edge #i contributes to H1
                    // (*edgeInfo)[i].second is already the zero column, so we change nothing
                    //(*edgeInfo)[i].second = H1pivotcounter;
                    H1pivotcounter++;
                    break;
                }
                pivot = std::max(nonpivot, otherNonpivot);
                nonpivot = std::min(nonpivot, otherNonpivot);
            }
            else{
                //haven't seen this pivot before, so store the reduced column and move on
                columns[pivot] = {nonpivot, pivot};
                output[pivot] = (*sortedDistances)[i].second;
                //the edge is dead (so can be ignored when computing H1)
                // we let (*edgeInfo)[i].first be a column with 1 in the row of pq and the rest 0
                //bit_tree_column columnWithOneAtpq;
                //columnWithOneAtpq.add_index(i);
                //i = (*edgeIndices)[p + q*(q - 1)/2]
                (*edgeInfo)[i].first = {i};
                (*edgeInfo)[i].second = (*sortedDistances)[i].first;
                //pivotsFound ++;
                break;
            }
        }
    }
    auto stopx = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_secondsx = stopx - startx;
    std::cout << "inside the H0 function: " << elapsed_secondsx.count() << " seconds\n";
    std::cout << "jCounter: " << jCounter << std::endl;
    //return output;
    return H1pivotcounter;
}

std::vector<std::pair<double, double>> H1barcode(const std::vector<std::vector<double>>* pointCloud){
    const int dim = static_cast<int>((*pointCloud)[0].size());
    const index noOfPoints = static_cast<index>((*pointCloud).size());
    const index noOfPairs = noOfPoints*(noOfPoints-1)/2;
    std::cout << noOfPairs << std::endl;
    column zeroColumn = {};
    std::vector<std::pair<column, double>> edgeInfo(noOfPairs, std::pair{zeroColumn, INFINITY});
    std::vector<std::vector<std::pair<double, index>>> distancesFromEachP;
    std::vector<std::pair<double, std::pair<index, index>>> sortedDistances;
    std::vector<index> edgeIndices(noOfPairs);
    const index noOfH1pivots = H0pivots(pointCloud, &edgeInfo, &distancesFromEachP, &sortedDistances, &edgeIndices);
    auto startx = std::chrono::steady_clock::now();
    std::vector<index> H1PivotEdges(noOfH1pivots);
    index pivotCounter = 0;
    std::vector<std::pair<double, double>> barcode(noOfH1pivots);
    return barcode;
    std::vector<bool> dead(noOfPairs, false);
    //return barcode;
    index counter = 0;
    index counterImportant = 0;
    index quickIterCounter = 0;
    index additionCounter1 = 0;
    index additionCounter2 = 0;
    index tripleCounter = 0;
    index zeroCounter = 0;
    index oblCounter = 0;
    index secondLoop = 0;
    index canSkipCounter = 0;
    bit_tree_column bitColumnSingleEdge, mxBitColumn, bitColumnpqr, tempCol, threshRedCol;
    bitColumnSingleEdge.init(noOfPairs);
    mxBitColumn.init(noOfPairs);
    bitColumnpqr.init(noOfPairs);
    for(index i=0; i < noOfPairs; i++){
        if(edgeInfo[i].first.size() == 3 && dead[edgeInfo[i].first[0]] && dead[edgeInfo[i].first[1]] && edgeInfo[i].second <= sortedDistances[i].first){
            tripleCounter ++;
            dead[i] = true;
            //bitColumnSingleEdge.add_index(i);
            //edgeInfo[i].first = {i};
            //edgeInfo[i].second = sortedDistances[i].first;
            continue;
            /* index altIndex = edgeInfo[i].first[0];
            if(edgeInfo[altIndex].first.size() == 1){
                continue;
            }
            //quickIterCounter ++;
            bitColumnSingleEdge.add_col(edgeInfo[altIndex].first);
            //shortcut = true; */
        }
        //bool shortcut = false;
        //bit_tree_column bitColumnSingleEdge;
        //bitColumnSingleEdge.init(noOfPairs);
        bitColumnSingleEdge.clear();
        bitColumnSingleEdge.add_index(i);
        //continue;
        //index altIndex = i;
        //reduce the edge to check if it's dead, or find its "death threshold"
        if(edgeInfo[i].first.size()){
            //bitColumn.add_col(edgeInfo[i].first);
            while(!bitColumnSingleEdge.is_empty()){
                const index mx = bitColumnSingleEdge.get_max_index();
                if(dead[mx]){
                    bitColumnSingleEdge.add_index(mx);
                    continue;
                }
                if(edgeInfo[mx].second > sortedDistances[i].first){
                    break;
                }
                if(edgeInfo[mx].first.size() > 1){
                    index mx2 = edgeInfo[mx].first[edgeInfo[mx].first.size()-2];
                    if(edgeInfo[mx2].second <= sortedDistances[i].first){
                        mxBitColumn.clear();
                        mxBitColumn.add_col(edgeInfo[mx].first);
                        mxBitColumn.add_col(edgeInfo[mx2].first);
                        mxBitColumn.get_col_and_clear(edgeInfo[mx].first);
                        if(edgeInfo[mx].first.size() == 1){
                            dead[mx] = true;
                        }
                    }
                }
                bitColumnSingleEdge.add_col(edgeInfo[mx].first);
            }
            if(bitColumnSingleEdge.is_empty()){
                dead[i] = true;
                /* edgeInfo[i].first = {i};
                //can I set this to 0?
                edgeInfo[i].second = sortedDistances[i].first; */
                quickIterCounter ++;
                /* if(shortcut){
                    edgeInfo[altIndex].first = {altIndex};
                } */
                continue;
            }
        }
        //else{
        //    bitColumnSingleEdge.add_index(i);
        //}
        // I think you can delete this without affecting correctness
        /* if(shortcut){
            bitColumnSingleEdge.get_col(edgeInfo[altIndex].first);
            edgeInfo[altIndex].first.push_back(altIndex);
        } */
        const index pIndex = sortedDistances[i].second.first;
        const index qIndex = sortedDistances[i].second.second;
        const std::vector<double> p = (*pointCloud)[pIndex];
        const std::vector<double> q = (*pointCloud)[qIndex];
        column columnSingleEdge = {};
        //index mxx = bitColumnSingleEdge.get_max_index();
        threshRedCol = bitColumnSingleEdge;
        index redmax = threshRedCol.get_max_index();
        bool importnt = bitColumnSingleEdge.get_max_index() == i;
        bitColumnSingleEdge.get_col_and_clear(columnSingleEdge);
        //go through the cofaces of edge i
        const double pq = sortedDistances[i].first;
        //------------------------------------
        //------------------------------------
        secondLoop ++;
        column acute(0);
        double threshold = INFINITY;
        for(const std::pair<double, index>& rPair : distancesFromEachP[pIndex]){
            const index rIndex = rPair.second;
            if(rIndex == qIndex){
                continue;
            }
            const std::vector<double> r = (*pointCloud)[rIndex];
            /* double dotProduct1 = 0;
            double dotProduct2 = 0;
            double dotProduct3 = 0;
            for (int d=0; d < dim; d++){
                dotProduct1 += (q[d] - p[d])*(r[d] - p[d]);
                dotProduct2 += (p[d] - q[d])*(r[d] - q[d]);
                dotProduct3 += (q[d] - r[d])*(p[d] - r[d]);
            }
            if(dotProduct1 <= 0 || dotProduct2 <= 0 || dotProduct3 <= 0){
                //the triangle has an angle >=90 deg, so not important
                oblCounter ++;
                continue;
            } */
           // This check seems pretty quick
            const int obt = obtuse(&r, &p, &q, dim);
            if(obt >= 2){
                //the triangle has an angle >=90 deg, so not important
                oblCounter ++;
                continue;
            }
            if(obt == 1){
                //the triangle has an angle >=90 deg, so not important
                /* if(acute.size() < 5){
                    acute.push_back(rIndex);
                }*/
                oblCounter ++;
                continue;
            }
            // Test for whether there is a point inside the circumsphere of p,q,r.
            // Has false negatives, but no false positives (I think)
            // Seems to skip only in a tiny fraction of cases, so I won't use it.
            /* bool canSkip = false;
            for (const index m : acute){
                const std::vector<double> s = (*pointCloud)[m];
                if(inCircumsphere3D(&p, &q, &r, &s)){
                    canSkip = true;
                    break;
                }
            }
            if(canSkip){
                continue;
            }
            if(acute.size() < 5){
                acute.push_back(rIndex);
            }*/
            canSkipCounter ++;
            if(importnt){
                counterImportant ++;
            }
            index smallpr = std::min(pIndex, rIndex);
            index largepr = std::max(pIndex, rIndex);
            index smallqr = std::min(qIndex, rIndex);
            index largeqr = std::max(qIndex, rIndex);
            const index prIndex = edgeIndices[smallpr + largepr*(largepr - 1)/2];
            const index qrIndex = edgeIndices[smallqr + largeqr*(largeqr - 1)/2];
            // pairwise distances between p, q, r
            const double pr = sortedDistances[prIndex].first;
            if(pr > threshold){
                break;
            }
            if(pr > 2*pq){
                counter ++;
            }
            while(redmax != -1 && pr > edgeInfo[redmax].second){
                threshRedCol.add_col(edgeInfo[redmax].first);
                redmax = threshRedCol.get_max_index();
            }
            if(redmax == -1){
                break;
            }
            const double qr = sortedDistances[qrIndex].first;
            // this is the square of the diameter of the circumcircle of pqr
            double pqrDiameter = 2*pq*pr*qr/(-(pq*pq+pr*pr+qr*qr)/2 + pq*pr+pr*qr+qr*pq);
            if(pqrDiameter > threshold){
                continue;
            }
            // initialize the column with nonzero at pq, pr, qr
            bitColumnpqr.clear();
            bitColumnpqr.add_index(prIndex);
            bitColumnpqr.add_index(qrIndex);
            while(!bitColumnpqr.is_empty()){
                const index mx = bitColumnpqr.get_max_index();
                if(dead[mx]){
                    additionCounter1 ++;
                    bitColumnpqr.add_index(mx);
                    continue;
                }
                if(!edgeInfo[mx].first.size() || edgeInfo[mx].second > pqrDiameter){
                    break;
                }
                bitColumnpqr.add_col(edgeInfo[mx].first);
            }
            if(bitColumnpqr.is_empty()){
                threshold = pqrDiameter;
            }
            bitColumnpqr.add_col(columnSingleEdge);
            while(!bitColumnpqr.is_empty()){
                const index mx = bitColumnpqr.get_max_index();
                if(dead[mx]){
                    additionCounter1 ++;
                    bitColumnpqr.add_index(mx);
                    continue;
                }
                if(!edgeInfo[mx].first.size()){
                    // we found a new pivot, so put the column into edgeInfo[mx]
                    column reducedColumn = {};
                    bitColumnpqr.get_col_and_clear(reducedColumn);
                    edgeInfo[mx].first = reducedColumn;
                    edgeInfo[mx].second = pqrDiameter;
                    H1PivotEdges[pivotCounter] = mx;
                    pivotCounter++;
                    break;
                }
                // < or <=?
                // add previous column with same pivot to this column
                //bitColumn.add_col(edgeInfo[mx].first);
                if(edgeInfo[mx].second <= pqrDiameter){
                    if(edgeInfo[mx].first.size() > 1){
                        index mx2 = edgeInfo[mx].first[edgeInfo[mx].first.size()-2];
                        if(edgeInfo[mx2].second <= sortedDistances[i].first){
                            mxBitColumn.clear();
                            mxBitColumn.add_col(edgeInfo[mx].first);
                            mxBitColumn.add_col(edgeInfo[mx2].first);
                            mxBitColumn.get_col_and_clear(edgeInfo[mx].first);
                            if(edgeInfo[mx].first.size() == 1){
                                dead[mx] = true;
                            }
                        }
                    }
                    additionCounter2 ++;
                    bitColumnpqr.add_col(edgeInfo[mx].first);
                }
                // a later column has the same pivot, so switch to reducing that one instead
                // if we want to keep track of pairing of birth and death simplices,
                // we also need to swap that here
                else{
                    const double temp = edgeInfo[mx].second;
                    edgeInfo[mx].second = pqrDiameter;
                    pqrDiameter = temp;
                    tempCol = bitColumnpqr;
                    bitColumnpqr.add_col(edgeInfo[mx].first);
                    tempCol.get_col_and_clear(edgeInfo[mx].first);
                }
            }
            if(bitColumnpqr.is_empty()){
                zeroCounter ++;
            }
        }
    }
    //std::cout << "Barcode with squared values:" << std::endl;
    for(index i=0; i < pivotCounter; i++){
        barcode[i] = {sortedDistances[H1PivotEdges[i]].first, edgeInfo[H1PivotEdges[i]].second};
        //std::cout << barcode[i].first << ", " << barcode[i].second << std::endl;
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
    std::cout << "counter:          " << counter << std::endl;
    std::cout << "canSkipCounter:   " << canSkipCounter << std::endl;
    std::cout << "counterImportant: " << counterImportant << std::endl;
    std::cout << "zeroCounter:      " << zeroCounter << std::endl;
    std::cout << "additionCounter1: " << additionCounter1 << std::endl;
    std::cout << "additionCounter2: " << additionCounter2 << std::endl;
    return barcode;
}


int main() {
    double min_coord = 0.0;
    double max_coord = 100.0;
    int dimensions = 3;
    int noOfPoints = 4000;
    
    std::vector<double> point = generateRandomPoint(min_coord, max_coord, dimensions);

    int noOfPairs = noOfPoints*(noOfPoints-1)/2;

    std::vector<std::vector<double>> pointCloud(noOfPoints, std::vector<double>(dimensions));
    for (int i=0; i < noOfPoints; i++){
        pointCloud[i] = generateRandomPoint(min_coord, max_coord, dimensions);
        //for (int j=0; j < dimensions; j++){
        //    std::cout << pointCloud[i][j] << ", ";
        //}
        //std::cout <<std::endl;
    }

    auto start = std::chrono::steady_clock::now();
    //std::vector<index> edgeIndices(noOfPairs);
    //std::vector<std::pair<double, std::pair<index, index>>> sD = allSortedDistances(&pointCloud, dimensions, &edgeIndices);
    //std::vector<std::pair<double, double>> barcode = H0barcode(&pointCloud, ...);
    //std::vector<std::pair<double, double>> barcode = 
    H1barcode(&pointCloud);
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
    auto mid1 = std::chrono::steady_clock::now();
    //vector<vector<int>> m = mu(&pointCloud, &asd);
    auto mid2 = std::chrono::steady_clock::now();
    //vector<vector<vector<pair<int, int>>>> phi = phiMatrix(&m);
    auto end = std::chrono::steady_clock::now();

    
    /*
    std::cout << std::endl;
    double avglength = 0;
    for (int j=0; j < numberOfPoints; j++){
        //std::cout << asd[j][0].second << std::endl;
        //std::cout << pointCloud[j][1] << std::endl;
        //std::cout << pointCloud[j][2] << std::endl << std::endl;
        for (int i=0; i < numberOfPoints; i++){
            avglength += phi[i][j].size();
        }
    }
    avglength = avglength/(numberOfPoints*numberOfPoints);
    cout << avglength << endl;
    */
    /*
     for (std::vector<double> p : pointCloud){
    //    for (std::vector<double> q : pointCloud){
    std::vector<std::vector<std::pair<double, int>>> allSortedDistances(numberOfPoints);
    for (int i=0; i < numberOfPoints; i++){
        //allSortedDistances.push_back(sortedDistances(pointCloud, i));
        allSortedDistances[i] = sortedDistances(&pointCloud, i);
        /* std::vector<std::pair<double, int>> srhiort = sortedDistances(pointCloud, i);
        for (int j=0; j < numberOfPoints; j++){
            //double distance = distpq(pointCloud[i], pointCloud[j]);
            double dx = pointCloud[i][0] - pointCloud[j][0];
            double dy = pointCloud[i][1] - pointCloud[j][1];
            double dz = pointCloud[i][2] - pointCloud[j][2];
            double distance = std::hypot(dx, dy, dz);
            //int x = 1;
        }
        
    }

    // Calculate the duration and cast to desired units (e.g., microseconds)
    std::chrono::duration<double> elapsed_seconds1 = mid1 - start;
    std::chrono::duration<double> elapsed_seconds2 = mid2 - mid1;
    std::chrono::duration<double> elapsed_seconds3 = end - mid2;

    std::cout << "First part: " << elapsed_seconds1.count() << " seconds\n";
    std::cout << "Second part: " << elapsed_seconds2.count() << " seconds\n";
    std::cout << "Second part: " << elapsed_seconds3.count() << " seconds\n";
    */
    return 0;
}