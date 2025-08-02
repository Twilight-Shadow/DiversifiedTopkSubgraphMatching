#include "ComputeSetIntersection.h"
#include <cstdint>

void ComputeSetIntersection::ComputeCandidates(const VertexID* larray, const ui l_count,
                                               const VertexID* rarray, const ui r_count,
                                               VertexID* cn, ui &cn_count) {
    if (l_count / 50 > r_count || r_count / 50 > l_count) {
        return ComputeCNGallopingAVX2(larray, l_count, rarray, r_count, cn, cn_count);
    }
    else {
        return ComputeCNMergeBasedAVX2(larray, l_count, rarray, r_count, cn, cn_count);
    }
}

void ComputeSetIntersection::ComputeCandidates(const VertexID* larray, const ui l_count,
                                               const VertexID* rarray, const ui r_count,
                                               ui &cn_count) {
        if (l_count / 50 > r_count || r_count / 50 > l_count) {
            return ComputeCNGallopingAVX2(larray, l_count, rarray, r_count, cn_count);
        }
        else {
            return ComputeCNMergeBasedAVX2(larray, l_count, rarray, r_count, cn_count);
        }
}

void ComputeSetIntersection::ComputeCNGallopingAVX2(const VertexID* larray, const ui l_count,
                                                    const VertexID* rarray, const ui r_count,
                                                    VertexID* cn, ui &cn_count) {
    cn_count = 0;

    if (l_count == 0 || r_count == 0)
        return;

    ui lc = l_count;
    ui rc = r_count;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        ui tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    ui li = 0;
    ui ri = 0;

    while (true) {
        while (larray[li] < rarray[ri]) {
            li += 1;
            if (li >= lc) {
                return;
            }
        }

        ri = GallopingSearchAVX2(rarray, ri, rc, larray[li]);
        if (ri >= rc) {
            return;
        }

        if (larray[li] == rarray[ri]) {
            cn[cn_count++] = larray[li];
            li += 1;
            ri += 1;
            if (li >= lc || ri >= rc) {
                return;
            }
        }
    }
}

void ComputeSetIntersection::ComputeCNGallopingAVX2(const VertexID* larray, const ui l_count,
                                                    const VertexID* rarray, const ui r_count,
                                                    ui &cn_count) {
    cn_count = 0;

    if (l_count == 0 || r_count == 0)
        return;

    ui lc = l_count;
    ui rc = r_count;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        ui tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    ui li = 0;
    ui ri = 0;

    while (true) {
        while (larray[li] < rarray[ri]) {
            li += 1;
            if (li >= lc) {
                return;
            }
        }

        ri = GallopingSearchAVX2(rarray, ri, rc, larray[li]);
        if (ri >= rc) {
            return;
        }

        if (larray[li] == rarray[ri]) {
            cn_count += 1;
            li += 1;
            ri += 1;
            if (li >= lc || ri >= rc) {
                return;
            }
        }
    }
}

void ComputeSetIntersection::ComputeCNMergeBasedAVX2(const VertexID* larray, const ui l_count,
                                                     const VertexID* rarray, const ui r_count,
                                                     VertexID* cn, ui &cn_count) {
    cn_count = 0;

    if (l_count == 0 || r_count == 0)
        return;

    ui lc = l_count;
    ui rc = r_count;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        ui tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    ui li = 0;
    ui ri = 0;

    __m256i per_u_order = _mm256_set_epi32(1, 1, 1, 1, 0, 0, 0, 0);
    __m256i per_v_order = _mm256_set_epi32(3, 2, 1, 0, 3, 2, 1, 0);
    VertexID* cur_back_ptr = cn;

    auto size_ratio = (rc) / (lc);
    if (size_ratio > 2) {
        if (li < lc && ri + 7 < rc) {
            __m256i u_elements = _mm256_set1_epi32(larray[li]);
            __m256i v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));

            while (true) {
                __m256i mask = _mm256_cmpeq_epi32(u_elements, v_elements);
                auto real_mask = _mm256_movemask_epi8(mask);
                if (real_mask != 0) {
                    // at most 1 element
                    *cur_back_ptr = larray[li];
                    cur_back_ptr += 1;
                }
                if (larray[li] > rarray[ri + 7]) {
                    ri += 8;
                    if (ri + 7 >= rc) {
                        break;
                    }
                    v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));
                } else {
                    li++;
                    if (li >= lc) {
                        break;
                    }
                    u_elements = _mm256_set1_epi32(larray[li]);
                }
            }
        }
    } else {
        if (li + 1 < lc && ri + 3 < rc) {
            __m256i u_elements = _mm256_loadu_si256((__m256i *) (larray + li));
            __m256i u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
            __m256i v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));
            __m256i v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);

            while (true) {
                __m256i mask = _mm256_cmpeq_epi32(u_elements_per, v_elements_per);
                auto real_mask = _mm256_movemask_epi8(mask);
                if (real_mask << 16 != 0) {
                    *cur_back_ptr = larray[li];
                    cur_back_ptr += 1;
                }
                if (real_mask >> 16 != 0) {
                    *cur_back_ptr = larray[li + 1];
                    cur_back_ptr += 1;
                }


                if (larray[li + 1] == rarray[ri + 3]) {
                    li += 2;
                    ri += 4;
                    if (li + 1 >= lc || ri + 3 >= rc) {
                        break;
                    }
                    u_elements = _mm256_loadu_si256((__m256i *) (larray + li));
                    u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
                    v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));
                    v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);
                } else if (larray[li + 1] > rarray[ri + 3]) {
                    ri += 4;
                    if (ri + 3 >= rc) {
                        break;
                    }
                    v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));
                    v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);
                } else {
                    li += 2;
                    if (li + 1 >= lc) {
                        break;
                    }
                    u_elements = _mm256_loadu_si256((__m256i *) (larray + li));
                    u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
                }
            }
        }
    }

    cn_count = (ui)(cur_back_ptr - cn);
    if (li < lc && ri < rc) {
        while (true) {
            while (larray[li] < rarray[ri]) {
                ++li;
                if (li >= lc) {
                    return;
                }
            }
            while (larray[li] > rarray[ri]) {
                ++ri;
                if (ri >= rc) {
                    return;
                }
            }
            if (larray[li] == rarray[ri]) {
                // write back
                cn[cn_count++] = larray[li];

                ++li;
                ++ri;
                if (li >= lc || ri >= rc) {
                    return;
                }
            }
        }
    }
    return;
}

void ComputeSetIntersection::ComputeCNMergeBasedAVX2(const VertexID* larray, const ui l_count,
                                                     const VertexID* rarray, const ui r_count,
                                                     ui &cn_count) {
    cn_count = 0;

    if (l_count == 0 || r_count == 0)
        return;

    ui lc = l_count;
    ui rc = r_count;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        ui tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    ui li = 0;
    ui ri = 0;

    constexpr int parallelism = 8;

    int cn_countv[parallelism] = {0, 0, 0, 0, 0, 0, 0, 0};
    __m256i sse_cn_countv = _mm256_load_si256((__m256i *) (cn_countv));
    __m256i sse_countplus = _mm256_set1_epi32(1);
    __m256i per_u_order = _mm256_set_epi32(1, 1, 1, 1, 0, 0, 0, 0);
    __m256i per_v_order = _mm256_set_epi32(3, 2, 1, 0, 3, 2, 1, 0);

    auto size_ratio = (rc) / (lc);
    if (size_ratio > 2) {
        if (li < lc && ri + 7 < rc) {
            __m256i u_elements = _mm256_set1_epi32(larray[li]);
            __m256i v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));

            while (true) {
                __m256i mask = _mm256_cmpeq_epi32(u_elements, v_elements);
                mask = _mm256_and_si256(sse_countplus, mask);
                sse_cn_countv = _mm256_add_epi32(sse_cn_countv, mask);
                if (larray[li] > rarray[ri + 7]) {
                    ri += 8;
                    if (ri + 7 >= rc) {
                        break;
                    }
                    v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));
                } else {
                    li++;
                    if (li >= lc) {
                        break;
                    }
                    u_elements = _mm256_set1_epi32(larray[li]);
                }
            }
            _mm256_store_si256((__m256i *) cn_countv, sse_cn_countv);
            for (int cn_countvplus : cn_countv) { cn_count += cn_countvplus; }
        }
    } else {
        if (li + 1 < lc && ri + 3 < rc) {
            __m256i u_elements = _mm256_loadu_si256((__m256i *) (larray + li));
            __m256i u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
            __m256i v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));
            __m256i v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);

            while (true) {
                __m256i mask = _mm256_cmpeq_epi32(u_elements_per, v_elements_per);
                mask = _mm256_and_si256(sse_countplus, mask);
                sse_cn_countv = _mm256_add_epi32(sse_cn_countv, mask);

                if (larray[li + 1] == rarray[ri + 3]) {
                    li += 2;
                    ri += 4;
                    if (li + 1 >= lc || ri + 3 >= rc) {
                        break;
                    }
                    u_elements = _mm256_loadu_si256((__m256i *) (larray + li));
                    u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
                    v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));
                    v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);
                } else if (larray[li + 1] > rarray[ri + 3]) {
                    ri += 4;
                    if (ri + 3 >= rc) {
                        break;
                    }
                    v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));
                    v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);
                } else {
                    li += 2;
                    if (li + 1 >= lc) {
                        break;
                    }
                    u_elements = _mm256_loadu_si256((__m256i *) (larray + li));
                    u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
                }
            }
        }
        _mm256_store_si256((__m256i *) cn_countv, sse_cn_countv);
        for (int cn_countvplus : cn_countv) { cn_count += cn_countvplus; }
    }

    if (li < lc && ri < rc) {
        while (true) {
            while (larray[li] < rarray[ri]) {
                ++li;
                if (li >= lc) {
                    return;
                }
            }
            while (larray[li] > rarray[ri]) {
                ++ri;
                if (ri >= rc) {
                    return;
                }
            }
            if (larray[li] == rarray[ri]) {
                cn_count++;
                ++li;
                ++ri;
                if (li >= lc || ri >= rc) {
                    return;
                }
            }
        }
    }
    return;
}

const ui ComputeSetIntersection::BinarySearchForGallopingSearchAVX2(const VertexID* array, ui offset_beg, ui offset_end, ui val) {
    while (offset_end - offset_beg >= 16) {
        auto mid = static_cast<uint32_t>((static_cast<unsigned long>(offset_beg) + offset_end) / 2);
        _mm_prefetch((char *) &array[(static_cast<unsigned long>(mid + 1) + offset_end) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &array[(static_cast<unsigned long>(offset_beg) + mid) / 2], _MM_HINT_T0);
        if (array[mid] == val) {
            return mid;
        } else if (array[mid] < val) {
            offset_beg = mid + 1;
        } else {
            offset_end = mid;
        }
    }

    // linear search fallback, be careful with operator>> and operation+ priority
    __m256i pivot_element = _mm256_set1_epi32(val);
    for (; offset_beg + 7 < offset_end; offset_beg += 8) {
        __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
        __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
        int mask = _mm256_movemask_epi8(cmp_res);
        if (mask != 0xffffffff) {
            return offset_beg + (_popcnt32(mask) >> 2);
        }
    }
    if (offset_beg < offset_end) {
        auto left_size = offset_end - offset_beg;
        __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
        __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
        int mask = _mm256_movemask_epi8(cmp_res);
        int cmp_mask = 0xffffffff >> ((8 - left_size) << 2);
        mask &= cmp_mask;
        if (mask != cmp_mask) { return offset_beg + (_popcnt32(mask) >> 2); }
    }
    return offset_end;
}

const ui ComputeSetIntersection::GallopingSearchAVX2(const VertexID* array, ui offset_beg, ui offset_end, ui val) {
    if (array[offset_end - 1] < val) {
        return offset_end;
    }

    // linear search
    __m256i pivot_element = _mm256_set1_epi32(val);
    if (offset_end - offset_beg >= 8) {
        __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
        __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
        int mask = _mm256_movemask_epi8(cmp_res);
        if (mask != 0xffffffff) { return offset_beg + (_popcnt32(mask) >> 2); }
    } else {
        auto left_size = offset_end - offset_beg;
        __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
        __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
        int mask = _mm256_movemask_epi8(cmp_res);
        int cmp_mask = 0xffffffff >> ((8 - left_size) << 2);
        mask &= cmp_mask;
        if (mask != cmp_mask) { return offset_beg + (_popcnt32(mask) >> 2); }
    }

    // galloping, should add pre-fetch later
    auto jump_idx = 8u;
    while (true) {
        auto peek_idx = offset_beg + jump_idx;
        if (peek_idx >= offset_end) {
            return BinarySearchForGallopingSearchAVX2(array, (jump_idx >> 1) + offset_beg + 1, offset_end, val);
        }
        if (array[peek_idx] < val) {
            jump_idx <<= 1;
        } else {
            return array[peek_idx] == val ? peek_idx :
                   BinarySearchForGallopingSearchAVX2(array, (jump_idx >> 1) + offset_beg + 1, peek_idx + 1, val);
        }
    }
}
