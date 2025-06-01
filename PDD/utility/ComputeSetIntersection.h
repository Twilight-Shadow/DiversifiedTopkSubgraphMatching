#ifndef COMPUTESETINTERSECTION_H
#define COMPUTESETINTERSECTION_H

#include "configuration/Types.h"
#include <immintrin.h>
#include <x86intrin.h>

/*
 * Because the set intersection is designed for computing common neighbors, the target is uieger.
 */

class ComputeSetIntersection {
public:
    static size_t galloping_cnt_;
    static size_t merge_cnt_;

    static void ComputeCandidates(const VertexID* larray, ui l_count, const VertexID* rarray,
                                  ui r_count, VertexID* cn, ui &cn_count);
    static void ComputeCandidates(const VertexID* larray, ui l_count, const VertexID* rarray,
                                  ui r_count, ui &cn_count);

    static void ComputeCNGallopingAVX2(const VertexID* larray, ui l_count,
                                       const VertexID* rarray, ui r_count, VertexID* cn,
                                       ui &cn_count);
    static void ComputeCNGallopingAVX2(const VertexID* larray, ui l_count,
                                       const VertexID* rarray, ui r_count, ui &cn_count);

    static void ComputeCNMergeBasedAVX2(const VertexID* larray, ui l_count, const VertexID* rarray,
                                        ui r_count, VertexID* cn, ui &cn_count);
    static void ComputeCNMergeBasedAVX2(const VertexID* larray, ui l_count, const VertexID* rarray,
                                        ui r_count, ui &cn_count);
    static const ui BinarySearchForGallopingSearchAVX2(const VertexID*  array, ui offset_beg, ui offset_end, ui val);
    static const ui GallopingSearchAVX2(const VertexID*  array, ui offset_beg, ui offset_end, ui val);
};


#endif // COMPUTESETINTERSECTION_H
