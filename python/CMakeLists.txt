# Python extension module nemo_astro._core (built with NEMO_PYTHON=ON,
# normally by pip through scikit-build-core; see pyproject.toml).

find_package(
  Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module
)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core MODULE src/bindings.cpp)
target_link_libraries(_core PRIVATE nemo)
target_compile_definitions(_core PRIVATE NEMO_PY_VERSION="${PROJECT_VERSION}")
target_compile_features(_core PRIVATE cxx_std_14)

install(TARGETS _core LIBRARY DESTINATION nemo_astro)
