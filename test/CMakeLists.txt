

#
# Add SandBox Testing 
# (Not included in Travis Testing)
#
add_subdirectory(sandbox/cdg/data)
add_subdirectory(sandbox/cdg/src)


#
# Add Travis Testing
#
add_subdirectory(travis/cdg)
add_subdirectory(travis/tbox)
