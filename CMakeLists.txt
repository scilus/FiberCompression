CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

cmake_policy(SET CMP0003 OLD)

PROJECT( FiberCompression )

#==================================================================
# Check if we are compiling on the mammoth
#==================================================================
OPTION(BUILD_FOR_MAMMOTH "Compile with paths specific to the Mammoth" OFF)
MARK_AS_ADVANCED( BUILD_FOR_MAMMOTH )

#==================================================================
# Option to build heavyweight version
#==================================================================
OPTION(BUILD_LIGHTWEIGHT "Build lightweight version" ON)
IF(BUILD_LIGHTWEIGHT)
	ADD_DEFINITIONS( -D_USE_LIGHTWEIGHT=1 )
ELSE(BUILD_LIGHTWEIGHT)
	ADD_DEFINITIONS( -D_USE_LIGHTWEIGHT=0 )
ENDIF(BUILD_LIGHTWEIGHT)
MARK_AS_ADVANCED(BUILD_LIGHTWEIGHT)

#==================================================================
# Option to save stats during process
#==================================================================
OPTION(BUILD_TO_SAVE_STATS "Build to save output statistics" OFF)
IF(BUILD_TO_SAVE_STATS)
	ADD_DEFINITIONS( -D_SAVE_STATS=1 )
ELSE(BUILD_TO_SAVE_STATS)
	ADD_DEFINITIONS( -D_SAVE_STATS=0 )
ENDIF(BUILD_TO_SAVE_STATS)
MARK_AS_ADVANCED(BUILD_TO_SAVE_STATS)

#==================================================================
# Look for QT4 (required)
#==================================================================
find_package(Qt4 4.6.0 REQUIRED)
set(QT_USE_QTGUI TRUE)
#set(QT_USE_QTOPENGL TRUE)
set(QT_USE_QTXML TRUE)
set(QT_USE_QTSVG TRUE)
include(${QT_USE_FILE})
include_directories(${QT_INCLUDE_DIR})

# LIBS
include_directories("libs")

# FFTW Configuration
IF(APPLE)
	find_library(FFTW NAMES fftw3 PATHS ${CMAKE_SOURCE_DIR}/libs/FFTW/OSX)
	set(FFTW_LIBRARIES "${FFTW}") 
ELSE(APPLE)
	IF(WIN32)
		set( FFTW_DIR "${CMAKE_SOURCE_DIR}/libs/FFTW/WIN32" CACHE FILEPATH FFTW_DIR )
		include_directories(${FFTW_DIR})
		link_directories(${FFTW_DIR})
	ELSE(WIN32)
		find_library(FFTW NAMES fftw3 PATHS ${CMAKE_SOURCE_DIR}/libs/FFTW/OSX)
		set(FFTW_LIBRARIES "${FFTW}") 
	ENDIF(WIN32)
ENDIF(APPLE)

# Nifti configuration
# Add this define for every platform. Used by nifti_io to try to read .nii.gz files.
ADD_DEFINITIONS(
    -DHAVE_ZLIB
)
file( GLOB NIFTI_HDR "${CMAKE_SOURCE_DIR}/libs/nifti/*.h" )
file( GLOB NIFTI_SRC "${CMAKE_SOURCE_DIR}/libs/nifti/*.c" )

# tclap lib
file( GLOB TCLAP "${CMAKE_SOURCE_DIR}/libs/tclap/*.h" )

# define sources
file( GLOB FC_SRC "src/*.cpp" )
list( REMOVE_ITEM FC_SRC "${CMAKE_SOURCE_DIR}/src/main_gui.cpp" "${CMAKE_SOURCE_DIR}/src/MainWindow.cpp" )

file( GLOB FC_HDR "src/*.h" )
list( REMOVE_ITEM FC_HDR "${CMAKE_SOURCE_DIR}/src/MainWindow.h" )

file( GLOB FC_SRC_GUI "src/*.cpp" )
list( REMOVE_ITEM FC_SRC_GUI "${CMAKE_SOURCE_DIR}/src/main.cpp" "${CMAKE_SOURCE_DIR}/src/Application.cpp" )
IF(WIN32)
	list( APPEND FC_SRC_GUI "${CMAKE_SOURCE_DIR}/rsc/images.rc" )
ENDIF(WIN32)

file( GLOB FC_HDR_GUI "src/*.h" )
list( REMOVE_ITEM FC_HDR_GUI "${CMAKE_SOURCE_DIR}/src/Application.h" )

file( GLOB FC_HDR_GUI_MOC "${CMAKE_SOURCE_DIR}/src/MainWindow.h" "${CMAKE_SOURCE_DIR}/src/Fibers.h" )
file( GLOB FC_HDR_MOC "${CMAKE_SOURCE_DIR}/src/Application.h" "${CMAKE_SOURCE_DIR}/src/Fibers.h" )
file( GLOB FC_RCC "${CMAKE_SOURCE_DIR}/rsc/images.qrc" )
set( exec_name fibercompression )
set( exec_gui_name fibercompression_gui )

# moc Qt source file without a ui file
QT4_WRAP_CPP( FC_MOC ${FC_HDR_MOC} ) 
QT4_WRAP_CPP( FC_GUI_MOC ${FC_HDR_GUI_MOC} ) 

# make rules to create c++ source from Qt resource files
QT4_ADD_RESOURCES(FC_RCC_SRC ${FC_RCC})

# Define include directories for other Igns libs
#include_directories( ${FC_BINARY_DIR} )     # needed to include qt-generated source files

#================================
# Define executables
#================================
SET( SOURCE ${FC_HDR} ${FC_SRC} ${FC_HDR_MOC} ${FC_MOC} ${NIFTI_SRC} ${NIFTI_HDR} ${TCLAP} )
ADD_EXECUTABLE( ${exec_name} ${executable_flags} ${SOURCE} )

SET( SOURCE_GUI ${FC_HDR_GUI} ${FC_HDR_GUI_MOC} ${FC_SRC_GUI} ${FC_GUI_MOC} ${FC_RCC_SRC} ${NIFTI_SRC} ${NIFTI_HDR} )
IF(APPLE)
	ADD_EXECUTABLE( ${exec_gui_name} MACOSX_BUNDLE ${executable_flags} ${SOURCE_GUI} )
	set( MACOSX_BUNDLE_ICON_FILE "${CMAKE_SOURCE_DIR}/rsc/icon.icns" )
	set( MACOSX_BUNDLE_BUNDLE_NAME "FiberCompression GUI" ) 
ELSE(APPLE)
	ADD_EXECUTABLE( ${exec_gui_name} ${executable_flags} ${SOURCE_GUI} )
ENDIF(APPLE)


TARGET_LINK_LIBRARIES( ${exec_name} ${QT_LIBRARIES} )
TARGET_LINK_LIBRARIES( ${exec_gui_name} ${QT_LIBRARIES} )

IF(FFTW_LIBRARIES)
	TARGET_LINK_LIBRARIES( ${exec_name} ${FFTW_LIBRARIES} )
	TARGET_LINK_LIBRARIES( ${exec_gui_name} ${FFTW_LIBRARIES} )
ELSE(FFTW_LIBRARIES)
	TARGET_LINK_LIBRARIES( ${exec_name} ${FFTW_DIR}/libfftw3-3.lib )
	TARGET_LINK_LIBRARIES( ${exec_gui_name} ${FFTW_DIR}/libfftw3-3.lib )
ENDIF(FFTW_LIBRARIES)

IF (APPLE)  
   TARGET_LINK_LIBRARIES( ${exec_name} z)
   TARGET_LINK_LIBRARIES( ${exec_gui_name} z)
ELSE(APPLE)
  IF(WIN32)
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES( ${exec_name} libz.so ) 
    TARGET_LINK_LIBRARIES( ${exec_gui_name} libz.so ) 
  ENDIF(WIN32)
ENDIF(APPLE)

# Need to specify the path for Mammoth, to overcome a bug in their setup.
IF (BUILD_FOR_MAMMOTH)
	set_target_properties( ${exec_name} PROPERTIES LINK_FLAGS "-L/opt/qt64/4.8.2/lib/")
ENDIF (BUILD_FOR_MAMMOTH)


MACRO(GetQtDLLs DEBUG_NAME RELEASE_NAME)
	FOREACH(module QT3SUPPORT QTOPENGL QTASSISTANT QTDESIGNER QTMOTIF QTNSPLUGIN
                   QAXSERVER QAXCONTAINER QTDECLARATIVE QTSCRIPT QTSVG
		QTUITOOLS QTHELP QTWEBKIT PHONON QTSCRIPTTOOLS QTMULTIMEDIA QTGUI QTTEST QTDBUS QTXML QTSQL QTXMLPATTERNS QTNETWORK QTCORE)
		IF(QT_USE_${module} OR QT_USE_${module}_DEPENDS)
			string(REPLACE "/lib/" "/bin/" QT_${module}_DLL "${QT_${module}_LIBRARY_DEBUG}")
			string(REPLACE ".lib" ".dll" QT_${module}_DLL "${QT_${module}_DLL}")
            LIST(APPEND ${DEBUG_NAME} "${QT_${module}_DLL}")
			string(REPLACE "/lib/" "/bin/" QT_${module}_DLL "${QT_${module}_LIBRARY_RELEASE}")
			string(REPLACE ".lib" ".dll" QT_${module}_DLL "${QT_${module}_DLL}")
            LIST(APPEND ${RELEASE_NAME} "${QT_${module}_DLL}")
        ENDIF()
	ENDFOREACH(module)
ENDMACRO()

IF(WIN32)
	# Copy qt dlls
    GetQtDLLs(DEBUG_DLLS RELEASE_DLLS)
    set (DLL_TO_DBG ${CMAKE_CURRENT_BINARY_DIR}/Debug)
    set (DLL_TO_RELEASE ${CMAKE_CURRENT_BINARY_DIR}/Release)

    foreach(dll ${DEBUG_DLLS})
        file(COPY ${dll} DESTINATION ${DLL_TO_DBG})
    endforeach()

    foreach(dll ${RELEASE_DLLS})
        file(COPY ${dll} DESTINATION ${DLL_TO_RELEASE})
    endforeach()

    install(FILES ${RELEASE_DLLS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} CONFIGURATIONS Release)
	
	# Copy FFTW dlls
	file( GLOB FFTW_DLLS "${FFTW_DIR}/*.dll" )
	foreach(dll ${FFTW_DLLS})
		file(COPY ${dll} DESTINATION ${DLL_TO_DBG})
        file(COPY ${dll} DESTINATION ${DLL_TO_RELEASE})
    endforeach()
	
	install(FILES ${FFTW_DLLS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} CONFIGURATIONS Release)
	
	# Copy .rc and .ico files
	file( GLOB RES "${CMAKE_SOURCE_DIR}/rsc/*.ico"
				   "${CMAKE_SOURCE_DIR}/rsc/*.png" 
				   "${CMAKE_SOURCE_DIR}/rsc/*.rc" )
	foreach(res ${RES})
		file(COPY ${res} DESTINATION ${DLL_TO_DBG})
        file(COPY ${res} DESTINATION ${DLL_TO_RELEASE})
    endforeach()
ENDIF(WIN32)
