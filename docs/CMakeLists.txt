# Set Doxygen input and output files.
set(DOXYGEN_INPUT_DIR
  ${PROJECT_SOURCE_DIR}/include
)
set(DOXYGEN_OUTPUT_DIR
  ${CMAKE_CURRENT_BINARY_DIR}/doxygen
)
set(DOXYFILE_IN
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
)
set(DOXYFILE_OUT
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
)

# Generate DoxyFile from the input file.
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

# Command for generating doc from Doxygen config file.
add_custom_target(docs
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating API documentation with Doxygen"
  VERBATIM
)
