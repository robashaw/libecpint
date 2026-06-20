add_executable(ref_atom_test test.cpp)
target_link_libraries(ref_atom_test ecpint)

add_test(NAME RefAtom
	 COMMAND ref_atom_test
	 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/lib/ref_atom)
