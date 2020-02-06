include_directories(${Hubbard_SOURCE_DIR}/include)

add_executable(hubbard-example Hubbard.cpp)
add_executable(anderson-example Anderson.cpp)
add_executable(anderson MyAnderson.cpp)
add_executable(anderson_chi MyAndersonChiSpin.cpp)
add_executable(anderson_chi_mats MyAndersonChiSpinMats.cpp)
add_executable(holstein-example HolsteinAnderson.cpp)

target_link_libraries(hubbard-example common-lib ${extlibs} ${parlibs})

target_link_libraries(anderson-example common-lib ${extlibs} ${parlibs})

target_link_libraries(anderson common-lib ${extlibs} ${parlibs})

target_link_libraries(anderson_chi common-lib ${extlibs} ${parlibs})

target_link_libraries(anderson_chi_mats common-lib ${extlibs} ${parlibs})

target_link_libraries(holstein-example common-lib ${extlibs} ${parlibs})
