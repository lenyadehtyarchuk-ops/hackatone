#----------------------------------------------------------------
# Generated CMake target import file for configuration "None".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "GDAL::GDAL" for configuration "None"
set_property(TARGET GDAL::GDAL APPEND PROPERTY IMPORTED_CONFIGURATIONS NONE)
set_target_properties(GDAL::GDAL PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_NONE "json-c::json-c;zstd::libzstd_shared;GEOS::geos_c;PROJ::proj;expat::expat;netCDF::netcdf"
  IMPORTED_LOCATION_NONE "${_IMPORT_PREFIX}/lib/x86_64-linux-gnu/libgdal.so.38.3.12.2"
  IMPORTED_SONAME_NONE "libgdal.so.38"
  )

list(APPEND _cmake_import_check_targets GDAL::GDAL )
list(APPEND _cmake_import_check_files_for_GDAL::GDAL "${_IMPORT_PREFIX}/lib/x86_64-linux-gnu/libgdal.so.38.3.12.2" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
