import distutils.core
import Cython.Build
distutils.core.setup(
    ext_modules = Cython.Build.cythonize("/media/paulosschlogl/Paulo/pauloscchlogl/Genome_kmers/kmer_counter/count_kmers.pyx"))
