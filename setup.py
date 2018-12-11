#############################################################################
#
#   Copyright 2018 Iain Mathieson
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#
#############################################################################

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext as build_pyx
import numpy

setup(name = "c_viterbi3",
      ext_modules=[Extension(name='c_viterbi3', sources=['c_viterbi3.pyx'], include_dirs=[numpy.get_include()])], cmdclass = {'build_ext': build_pyx})

