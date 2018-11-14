

include( ExternalProject )

set( boost_URL "https://dl.bintray.com/boostorg/release/1.68.0/source/boost_1_68_0.tar.gz" )
set( boost_SHA256 "da3411ea45622579d419bfda66f45cd0f8c32a181d84adfa936f5688388995cf" )
set( boost_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/third_party/boost )
set( boost_BUILD   ${boost_INSTALL}/build )
set( boost_INCLUDE_DIR ${boost_BUILD}/include )
set( boost_LIB_DIR ${boost_BUILD}/lib )

ExternalProject_Add( boost
        PREFIX boost
        URL ${boost_URL}
        URL_HASH SHA256=${boost_SHA256}
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND
        ./bootstrap.sh
        --with-libraries=filesystem
        --with-libraries=mpi
        --with-libraries=serialization
        --prefix=${boost_BUILD}
        BUILD_COMMAND
        ./b2 install link=static variant=release threading=multi runtime-link=static
        INSTALL_COMMAND ""
        INSTALL_DIR ${boost_INSTALL} )

set( Boost_LIBRARIES
        ${boost_LIB_DIR}/libboost_filesystem.a
        ${boost_LIB_DIR}/libboost_mpi.a
        ${boost_LIB_DIR}/libboost_serialization.a )
message( STATUS "Boost static libs: " ${Boost_LIBRARIES})

