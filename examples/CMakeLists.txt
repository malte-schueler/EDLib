
add_executable(hubbard-example Hubbard.cpp)
add_executable(anderson_imag Anderson.cpp)
add_executable(anderson_real Anderson_real.cpp )
add_executable(anderson_chi_imag Anderson_chi_mats.cpp)
add_executable(holstein-example HolsteinAnderson.cpp)

target_link_libraries(hubbard-example common-lib EDLib ${extlibs} ${parlibs})

target_link_libraries(anderson_imag common-lib EDLib ${extlibs} ${parlibs})

target_link_libraries(anderson_real common-lib EDLib ${extlibs} ${parlibs})

target_link_libraries(anderson_chi_imag common-lib EDLib ${extlibs} ${parlibs})

target_link_libraries(holstein-example common-lib EDLib ${extlibs} ${parlibs})
