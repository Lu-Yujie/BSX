#ifndef SETOP_H
#define SETOP_H
/**
 * some operations on set
 * 24.2.23
*/
// multi sets union
#include <iostream>
#include <vector>
#include <queue>
#include <climits>
using namespace std;

typedef uint32_t ui;

namespace SetOp {

struct ArrayElement {
    ui value;
    ui arrayIndex;
    ui elementIndex;

    ArrayElement(ui v, ui array_idx, ui ele_idx): value(v), arrayIndex(array_idx), elementIndex(ele_idx) {}

    // min-heap
    bool operator>(const ArrayElement& other) const {
        return value > other.value;
    }
};

vector<ui> unionMultiple(const vector<vector<ui>>& arrays) {
    vector<ui> result;
    priority_queue<ArrayElement, vector<ArrayElement>, greater<ArrayElement>> pq;

    // init-pq
    for (ui i = 0; i < arrays.size(); ++i) {
        if (!arrays[i].empty()) {
            pq.push({arrays[i][0], i, 0});
        }
    }

    if (pq.empty()) return result;

    // first element
    ArrayElement current = pq.top();
    pq.pop();
    result.push_back(current.value);
    if (current.elementIndex + 1 < arrays[current.arrayIndex].size()) {
        pq.push({arrays[current.arrayIndex][current.elementIndex + 1], current.arrayIndex, current.elementIndex + 1});
    }

    // scan all array
    while (!pq.empty()) {
        ArrayElement current = pq.top();
        pq.pop();

        if (result.back() != current.value) {
            result.push_back(current.value);
        }

        // 如果当前数组还有元素，将下一个元素加入到优先队列
        if (current.elementIndex + 1 < arrays[current.arrayIndex].size()
            // a possible opt, if next value in result.back, push the next one (if rarely happens, it's neg-opt)
            // && result.back() != arrays[current.arrayIndex][current.elementIndex + 1]
            ) {
            pq.push({arrays[current.arrayIndex][current.elementIndex + 1], current.arrayIndex, current.elementIndex + 1});
        }
    }

    return result;
}

vector<ui> unionMultiple(const vector<vector<ui>*>& arrays) {
    vector<ui> result;
    priority_queue<ArrayElement, vector<ArrayElement>, greater<ArrayElement>> pq;

    // init-pq
    for (ui i = 0; i < arrays.size(); ++i) {
        if (!(*(arrays[i])).empty()) {
            pq.push({(*(arrays[i]))[0], i, 0});
        }
    }

    if (pq.empty()) return result;

    // first element
    ArrayElement current = pq.top();
    pq.pop();
    result.push_back(current.value);
    if (current.elementIndex + 1 < (*(arrays[current.arrayIndex])).size()) {
        pq.push({(*(arrays[current.arrayIndex]))[current.elementIndex + 1], current.arrayIndex, current.elementIndex + 1});
    }

    // scan all array
    while (!pq.empty()) {
        ArrayElement current = pq.top();
        pq.pop();

        if (result.back() != current.value) {
            result.push_back(current.value);
        }

        // 如果当前数组还有元素，将下一个元素加入到优先队列
        if (current.elementIndex + 1 < (*(arrays[current.arrayIndex])).size()) {
            pq.push({(*(arrays[current.arrayIndex]))[current.elementIndex + 1], current.arrayIndex, current.elementIndex + 1});
        }
    }

    return result;
}

vector<ui> unionMultiple(const ui** arrays, const ui* arrays_size, const ui arrays_num) {
    vector<ui> result;
    priority_queue<ArrayElement, vector<ArrayElement>, greater<ArrayElement>> pq;

    // init-pq
    for (ui i = 0; i < arrays_num; ++i) {
        if (arrays_size[i] != 0) {
            pq.push({arrays[i][0], i, 0});
        }
    }

    if (pq.empty()) return result;

    // first element
    ArrayElement current = pq.top();
    pq.pop();
    result.push_back(current.value);
    if (current.elementIndex + 1 < arrays_size[current.arrayIndex]) {
        pq.push({arrays[current.arrayIndex][current.elementIndex + 1], current.arrayIndex, current.elementIndex + 1});
    }

    // scan all array
    while (!pq.empty()) {
        ArrayElement current = pq.top();
        pq.pop();

        if (result.back() != current.value) {
            result.push_back(current.value);
        }

        // 如果当前数组还有元素，将下一个元素加入到优先队列
        if (current.elementIndex + 1 < arrays_size[current.arrayIndex]
            // a possible opt, if next value in result.back, push the next one (if rarely happens, it's neg-opt)
            // && result.back() != arrays[current.arrayIndex][current.elementIndex + 1]
            ) {
            pq.push({arrays[current.arrayIndex][current.elementIndex + 1], current.arrayIndex, current.elementIndex + 1});
        }
    }

    // ui* returned = new ui[result.size()];
    // copy(result.begin(), result.end(), returned);
    return result;
}

void unionTwoAndUpdate(vector<ui>& array1, vector<ui>& array2) {
    vector<ui> unions;
    ui i = 0, j = 0;
    ui array1_size = array1.size();
    ui array2_size = array2.size();
    while (i < array1_size && j < array2_size) {
        if (array1[i] < array2[j]) {
            unions.push_back(array1[i]);
            ++i;
        } else if (array1[i] > array2[j]) {
            unions.push_back(array2[j]);
            ++j;
        } else {
            unions.push_back(array1[i]);
            ++i;
            ++j;
        }
    }

    // Add remaining elements from array1
    while (i < array1_size) {
        unions.push_back(array1[i]);
        ++i;
    }
    // Add remaining elements from array2
    while (j < array2_size) {
        unions.push_back(array2[j]);
        ++j;
    }
    array1.swap(unions);
}

bool setContain(const vector<ui>& big, const vector<ui>& small) {
    ui i = 0, j = 0;
    ui big_size = big.size(), small_size = small.size();

    while (i < big_size && j < small_size) {
        if (big[i] == small[j]) {
            ++i; ++j;
        } else if (big[i] > small[j]) {
            return false;
        } else {
            ++i;
        }
    }

    return j == small_size;
}

bool setContain(const ui* big, ui big_size, const ui* small, ui small_size) {
    ui i = 0, j = 0;

    while (i < big_size && j < small_size) {
        if (big[i] == small[j]) {
            ++i; ++j;
        } else if (big[i] > small[j]) {
            return false;
        } else {
            ++i;
        }
    }

    return j == small_size;
}

vector<ui> intersectTwo(const ui* array1, const ui* array2, ui array1_size, ui array2_size) {
    vector<ui> intersection;
    ui i = 0, j = 0;
    while (i < array1_size && j < array2_size) {
        if (array1[i] < array2[j]) {
            ++i;
        } else if (array1[i] > array2[j]) {
            ++j;
        } else {
            intersection.push_back(array1[i]);
            ++i;
            ++j;
        }
    }
    return intersection;
}

vector<ui> intersectTwo(const vector<ui>& array1, const vector<ui>& array2) {
    vector<ui> intersection;
    ui i = 0, j = 0;
    while (i < array1.size() && j < array2.size()) {
        if (array1[i] < array2[j]) {
            ++i;
        } else if (array1[i] > array2[j]) {
            ++j;
        } else {
            intersection.push_back(array1[i]);
            ++i;
            ++j;
        }
    }
    return intersection;
}

vector<ui> intersectTwo(const vector<ui>& array1, const ui* array2, ui array2_size) {
    vector<ui> intersection;
    ui i = 0, j = 0;
    while (i < array1.size() && j < array2_size) {
        if (array1[i] < array2[j]) {
            ++i;
        } else if (array1[i] > array2[j]) {
            ++j;
        } else {
            intersection.push_back(array1[i]);
            ++i;
            ++j;
        }
    }
    return intersection;
}

bool haveOverlapTwo(const vector<ui>& array1, const vector<ui>& array2) {
    ui i = 0, j = 0;
    while (i < array1.size() && j < array2.size()) {
        if (array1[i] < array2[j]) {
            ++i;
        } else if (array1[i] > array2[j]) {
            ++j;
        } else {
            ++i;
            ++j;
            return true;
        }
    }
    return false;
}

vector<ui> intersectMultiple(const ui** arrays, const ui* arrays_size, const ui arrays_num) {
    vector<ui> result;

    // check arrays.size
    if (arrays_num == 0) return result;
    if (arrays_num == 1) {
        result.reserve(arrays_size[0]);
        for (ui i = 0; i < arrays_size[0]; i++) result.emplace_back(arrays[0][i]);
        return result;
    }

    // init pointer array
    vector<ui> pointers(arrays_num, 0);
    pointers[0] = 1;
    if (arrays_size[0] == 0) return result;
    ui minVal = arrays[0][0];
    ui cnt = 1;
    ui array_idx = 1;

    // until pointers[*] >= arrays[*].size(), aka, one array reach the end
    while (true) {
        if (pointers[array_idx] >= arrays_size[array_idx]) {
            return result;
        }
        while (arrays[array_idx][pointers[array_idx]] < minVal) {
            pointers[array_idx]++;
            if (pointers[array_idx] >= arrays_size[array_idx]) {
                return result;
            }
        }
        if (arrays[array_idx][pointers[array_idx]] != minVal) {
            minVal = arrays[array_idx][pointers[array_idx]];
            cnt = 1;
        } else {
            cnt++;
            if (cnt == arrays_num) {
                result.emplace_back(minVal);
                // next line is right, but whether adding it depends on the expectation of #result
                //   if expectation of #result is small, aka, this check is less than arrays_num, then it worths.
                //   else it does't worth
                // if (pointers[array_idx] == arrays_size[array_idx] - 1) return result;
            }
        }
        pointers[array_idx]++;
        array_idx++;
        array_idx%=arrays_num;
    }
    return result;
}

// A = A-B
void setDifference(vector<ui>& A,vector<ui>& B) {
    ui i = 0;  // idx_A
    ui j = 0;  // idx_B
    ui k = 0;  // idx_result, also size
    auto A_size = A.size();
    auto B_size = B.size();

    while (i < A_size && j < B_size) {
        if (A[i] < B[j]) {  // A[i] not in B, add to result(k)
            A[k++] = A[i++];
        }
        else if (A[i] > B[j]) {  // A[i] > B[j], next ele of B
            j++;
        }
        else {  // A[i] == B[j], delete A[i]
            i++;
        }
    }

    // add remained eles of A
    while (i < A_size) {
        A[k++] = A[i++];
    }
    A.resize(k);
}

void setDifference(ui* A, ui& A_size, vector<ui>&B) {
    ui i = 0;  // idx_A
    ui j = 0;  // idx_B
    ui k = 0;  // idx_result, also size
    ui B_size = B.size();

    while (i < A_size && j < B_size) {
        if (A[i] < B[j]) {  // A[i] not in B, add to result(k)
            A[k++] = A[i++];
        }
        else if (A[i] > B[j]) {  // A[i] > B[j], next ele of B
            j++;
        }
        else {  // A[i] == B[j], delete A[i]
            i++;
        }
    }

    // add remained eles of A
    while (i < A_size) {
        A[k++] = A[i++];
    }
    A_size = k;
}

// A = A∩B
ui intersectAndUpdate(vector<VertexID>&A, vector<VertexID>&B) {
    ui indexA = 0;
    ui indexB = 0;
    ui insertPos = 0;
    ui A_size = A.size();
    ui B_size = B.size();

    while (indexA < A_size && indexB < B_size) {
        if (A[indexA] < B[indexB]) {
            indexA++;
        }
        else if (A[indexA] > B[indexB]) {
            indexB++;
        }
        else {
            A[insertPos++] = A[indexA];
            indexA++;
            indexB++;
        }
    }

    A.resize(insertPos);
    return insertPos;
}


} // namespace SetOp

#endif
