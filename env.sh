export SCRAM_ARCH=slc6_amd64_gcc630
export CMSSW_RELEASE=CMSSW_9_3_9

function loadtool() {
  if ! [ $CMSSW_RELEASE ]
  then
    echo "No release set."
    return
  fi

  PACKAGE=$1

  CMSSW_RELEASE_BASE=/cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/cmssw/$CMSSW_RELEASE
  TOOLBOX=$CMSSW_RELEASE_BASE/config/toolbox/$SCRAM_ARCH/tools/selected

  [ -e ${TOOLBOX}/${PACKAGE}.xml ] || return

  echo source $(sed -n 's/.*environment name="[^"]*_BASE" *default="\([^"]*\)".*/\1/p' ${TOOLBOX}/${PACKAGE}.xml)/etc/profile.d/init.sh
  source $(sed -n 's/.*environment name="[^"]*_BASE" *default="\([^"]*\)".*/\1/p' ${TOOLBOX}/${PACKAGE}.xml)/etc/profile.d/init.sh
}

for TOOL in gcc-cxxcompiler glibc python pcre root_interface xz libtiff libjpg libpng boost curl openssl hepmc lhapdf fastjet tbb
do
  loadtool $TOOL
done
export LHAPATH=$(ls /cvmfs/cms.cern.ch/$SCRAM_ARCH/external/lhapdf/*/share/LHAPDF | head -n1)
