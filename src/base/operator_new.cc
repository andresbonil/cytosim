// Cytosim was created by Francois Nedelec. Copyright 2019 Cambridge University

// replace global new & delete operators to control memory alignment

#include <new>
#include <cstdio>
#include <cstdlib>

void* operator new(std::size_t size)
{
    //printf("new(%lu)\n", size);
    void * ptr = nullptr;
#if ( 1 )
    constexpr std::size_t sup = 1 << 30;
    // we align all memory to 32 bytes
    if ( size < sup )
    {
        if ( posix_memalign(&ptr, 32, size) )
            throw std::bad_alloc();
    }
#else
    // system's default (unaligned) memory
    ptr = std::malloc(size);
#endif
    if ( ptr == nullptr )
        throw std::bad_alloc();
    //std::printf("Cytosim new %5zu %p\n", s, ptr);
    return ptr;
}


void operator delete(void * ptr) throw()
{
    //std::printf("Cytosim delete    %p\n", ptr);
    std::free(ptr);
}


/*
void* operator new[](std::size_t s) throw(std::bad_alloc)
{
    std::printf("Cytosim new[] %5zu\n", s);
    return ::operator new(s);
}

void operator delete[](void *ptr) throw()
{
    std::printf("Cytosim delete[]    %p\n", ptr);
    ::operator delete(ptr);
}
*/
