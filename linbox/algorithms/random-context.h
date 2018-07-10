/* Copyright (C) 2018 the members of the LinBox group
 * Written by Daniel S. Roche <roche@usna.edu>
 *
 * This file is part of the LinBox library.
 *
 * ========LICENCE========
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * LinBox is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.     See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_RANDOM_CONTEXT_H
#define __LINBOX_RANDOM_CONTEXT_H

/** @file algorithms/random_context.h
 * @brief Randomization context for algorithm execution
 */

#include <utility>
#include <memory>

/* choices: Heuristic, MonteCarlo, LasVegas, Deterministic */
#ifndef _LINBOX_DEFAULT_RANDOM_CAT
#define _LINBOX_DEFAULT_RANDOM_CAT Heuristic
#endif

#define LB_MACRO_GLUE2(x,y) x ## y
#define LB_MACRO_GLUE(x,y) LB_MACRO_GLUE2(x,y)
#define LB_DEF_RC LB_MACRO_GLUE(_LINBOX_DEFAULT_RANDOM_CAT,Context)

namespace LinBox
{

    struct DeterministicAlgorithm {} ;

    struct LasVegasAlgorithm {
        using StricterAlgorithm = DeterministicAlgorithm;
        static constexpr StricterAlgorithm stricter()
        { return StricterAlgorithm(); }
    };

    struct MonteCarloAlgorithm {
        using StricterAlgorithm = LasVegasAlgorithm;
        static constexpr StricterAlgorithm stricter()
        { return StricterAlgorithm(); }
    };

    struct HeuristicAlgorithm {
        using StricterAlgorithm = MonteCarloAlgorithm;
        static constexpr StricterAlgorithm stricter()
        { return StricterAlgorithm(); }
    };


    template <typename RNG_Type>
    struct RNGContext {
        using RNG = RNG_Type;
        using Self = RNGContext<RNG>;

    protected :
        std::shared_ptr<RNG> _rng;

    public :
        RNGContext(std::shared_ptr<RNG> rng) :
            _rng(std::move(rng))
        {}

        template <typename... Args>
        RNGContext(Args&&... args) :
            _rng(new RNG(std::forward<Args>(args)...))
        {}

        // use the default copy and move constructors
        RNGContext(const Self&) = default;
        RNGContext(Self&&) = default;

        RNG& generator() { return *_rng; }

        decltype((*_rng)()) operator() ()
        {
            return (*_rng)();
        }
    };


    struct DeterministicContext
    {
        using Tag = DeterministicAlgorithm;
        static constexpr Tag getTag() { return Tag(); }
    };


    template <typename RNG_Type>
    struct LasVegasContextRNG :public RNGContext<RNG_Type> {
        using RNG = RNG_Type;
        using Parent = RNGContext<RNG>;
        using Tag = LasVegasAlgorithm;

        using Parent::RNGContext; // inherit RNGContext constructors

        static constexpr Tag getTag() { return Tag(); }
    };
    using LasVegasContext = LasVegasContextRNG<Givaro::GivRandom>;


    template <typename RNG_Type>
    struct MonteCarloContextRNG :public RNGContext<RNG_Type> {
        using RNG = RNG_Type;
        using Parent = RNGContext<RNG>;
        using Self = MonteCarloContextRNG<RNG>;
        using Tag = MonteCarloAlgorithm;
        static constexpr double DEFAULT_ERRBOUND = 0.000001;

    private :
        double _errbound = DEFAULT_ERRBOUND;
        int _weight = 0; // how many levels deep of a copy is this
        bool _original = true; // false if this was created by copy constructor

    public :
        MonteCarloContextRNG() {}

        template <typename... RNGArgs>
        MonteCarloContextRNG(double errbound, RNGArgs&&... rng_args) :
            Parent(std::forward<RNGArgs>(rng_args)...),
            _errbound(errbound)
        {}

        MonteCarloContextRNG(double errbound, std::shared_ptr<RNG> rng) :
            Parent(std::move(rng)),
            _errbound(errbound)
        {}

        /** @brief Construct a new Context with a split of the error probability.
         * We try to split in such a way that repeated splits don't decrease
         * the error probability exponentially (if the repeated splits
         * are all the same "direction" of copying).
         */
        MonteCarloContextRNG (Self& other) :
            Parent(other._rng),
            _original(false)
        {
            if (other._weight == 0) {
                // first copy; split evenly and both get weight 1
                _errbound = other._errbound /= 2;
                _weight = other._weight = 1;
            }
            else {
                // divide other._errbound into smaller and larger portions
                auto lightp = other._errbound / (other._weight + 2);
                auto heavyp = lightp * (other._weight + 1);
                if (other._original) {
                    // other has already been copied from; it keeps more probability
                    other._errbound = heavyp;
                    other._weight ++;
                    _errbound = lightp;
                    _weight = 0;
                }
                else {
                    // other was the copy; more probability to this copy
                    _errbound = heavyp;
                    _weight = other._weight + 1;
                    other._errbound = lightp;
                    other._weight = 0;
                }
            }
            other._original = true;
        }

        /** @brief Constructs a new Context taking ALL of this one's error probability.
         */
        MonteCarloContextRNG (Self&& other) :
            Parent(std::move(other._rng)),
            _errbound(other._errbound),
            _weight(other._weight),
            _original(other._original)
        {
            other._errbound = 0.;
            other._weight = -1;
        }

        // disable assignment and const copying
        MonteCarloContextRNG (const Self&) = delete;
        Self& operator= (const Self&) = delete;
        Self& operator= (Self&&) = delete;

        /* @brief Creates a new Context object with a fraction of this error probability.
         */
        Self&& split (double result_fraction = 1.) {
            Self res(_errbound * result_fraction, this->_rng);
            _errbound -= res._errbound;
            _weight = 0;
            return std::move(res);
        }

        constexpr double errbound() const { return _errbound; }

        static constexpr Tag getTag() { return Tag(); }
    };
    using MonteCarloContext = MonteCarloContextRNG<Givaro::GivRandom>;


    template <typename RNG_Type>
    struct HeuristicContextRNG :public RNGContext<RNG_Type> {
        using RNG = RNG_Type;
        using Parent = RNGContext<RNG>;
        using Tag = HeuristicAlgorithm;

        using Parent::RNGContext; // inherit RNGContext constructors

        static constexpr double errbound() { return 0.5; }

        static constexpr Tag getTag() { return Tag(); }
    };
    using HeuristicContext = HeuristicContextRNG<Givaro::GivRandom>;

#if 0
    // Here is a generic declaration of a function that takes one
    // actual argument of type T.
    // First there is a generic version that tries to move to
    // a stricter and stricter probability class until an implemention
    // is found.
    // The second thing is just a convenience function for the initial
    // call to the method that fills in default arguments.
    // This is boilerplate that has to be copied for every funtion.

    template <typename T, class RandomContext, class RandomTag>
    inline int foo(const T& data, RandomContext&& context, RandomTag tag)
    { return foo(data, std::forward<RandomContext>(context), tag.stricter()); }

    template <typename T, class RandomContext>
    inline int foo(const T& data, RandomContext&& context)
    { return foo(data, std::forward<RandomContext>(context), context.getTag); }

    // To provide an actual implementation, the context should be passed by value
    // and the tag should be specified for that implementation's type of algorithm.
    // Here are two implementations of foo, first a Las Vegas algorithm, and then
    // a Monte Carlo algorithm.

    template <typename T, class RandomContext>
    int foo(const T& data, RandomContext rc, LasVegasAlgorithm) {
        do { x = data + rc() }
        while (x == 10);
        return x;
    }

    template <typename T, class RandomContext>
    int foo(const T& data, RandomContext rc, MonteCarloAlgorithm) {
        return data + (rc() % (int)(1 / rc.errbound()));
    }
#endif

    LB_DEF_RC& defaultRandomContext() {
        static LB_DEF_RC rc;
        return rc;
    }

} //LinBox

#undef LB_DEF_RC

#endif //__LINBOX_RANDOM_CONTEXT_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
