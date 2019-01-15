#pragma once

namespace LinBox {

    /**
     * \brief Mono-denominator rational vector.
     *
     * Stores a rational vector of field element as
     * a vector of numerators and one scalar denominator.
     */
    template <class Vector>
    struct MdrVector {
        using Field = typename Vector::Field;
        using Element = typename Field::Element;

        // @fixme Do we need accesors to elements though methods,
        // so that it can be used like any vector?
        // If so, there won't be any reference to element.

        Vector num;
        Element den;

        MdrVector(const Field& F, size_t n)
            : num(F, n)
        {
        }

        size_t size() const { return num.size(); }
    };
}
