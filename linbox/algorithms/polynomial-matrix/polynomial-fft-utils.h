#ifndef __LINBOX_polynomial_fft_utils_H
#define __LINBOX_polynomial_fft_utils_H

#include <type_traits>

namespace LinBox {

    namespace FFT_utils {

        //template <class Field, class Simd>
        //using is_compatible_simd_v = std::is_same<typename Field::Element, typename Simd::scalar_t>::value;

        template< bool B, class T = void >
        using enable_if_t = typename std::enable_if<B, T>::type;

        template <class Field, class T = void>
        using enable_if_integral_field_t = enable_if_t<Field::is_elt_integral_v, T>;

        template <class Field, class T = void>
        using enable_if_floating_field_t = enable_if_t<Field::is_elt_floating_point_v, T>;

        template <class Field, class Simd, class T = void>
        using enable_if_integral_compatible_simd_t = enable_if_t<Field::is_elt_integral_v && Simd::template is_same_element<Field>::value, T>;

        template <class Field, class Simd, class T = void>
        using enable_if_floating_compatible_simd_t = enable_if_t<Field::is_elt_floating_point_v && Simd::template is_same_element<Field>::value, T>;

        template <class Field, class Simd, class T = void>
        using enable_if_same_element_t = enable_if_t<Simd::template is_same_element<Field>::value, T>;

        template <class T1, class T2, class T = void>
        using enable_if_same_t = enable_if_t<std::is_same<T1, T2>::value, T>;
    }
}
#endif // __LINBOX_polynomial_fft_utils_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
