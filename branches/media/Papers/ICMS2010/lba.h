#ifdef KAAPI_STL
namespace LinBox {
  namespace lba {
    using kstl::transform;
    using kstl::for_each;
    using kstl::copy;
    ...
  }
}
#elsif defined(OMP)
namespace LinBox {
  namespace lba {
    using _gnu_parallel::transform;
    using _gnu_parallel::for_each;
    using _gnu_parallel::copy;
    ...
  }
}
#else
namespace LinBox {
  namespace lba {
    using std::transform;
    using std::for_each;
    using std::copy;
    ...
  }
}

#endif
