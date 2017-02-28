FROM ubuntu:14.04

MAINTAINER Jason W. DeGraw jason.degraw@nrel.gov

# Github repo
ENV ENERGYPLUS_REPO_URL https://github.com/NREL/EnergyPlus.git


# Set up, clone the repo, and build
RUN apt-get update && apt-get install -y git cmake g++ gfortran python \
    && mkdir energyplus \
    && cd energyplus \
    && git clone -b develop --single-branch $ENERGYPLUS_REPO_URL . \
    && export ENERGYPLUS_SHA=$(git rev-parse --short=10 HEAD) \
    && export MAJOR=$(grep -Po 'set\( CMAKE_VERSION_MAJOR \K[0-9] (?=\))' CMakeLists.txt) \
    && export MINOR=$(grep -Po 'set\( CMAKE_VERSION_MINOR \K[0-9] (?=\))' CMakeLists.txt) \
    && export PATCH=$(grep -Po 'set\( CMAKE_VERSION_PATCH \K[0-9] (?=\))' CMakeLists.txt) \
    && export ENERGYPLUS_VERSION=$(echo $MAJOR.$MINOR.$PATCH) \
    && mkdir build \
    && cd build \
    && cmake -DBUILD_TESTING=ON -DBUILD_PACKAGE=ON -DCPACK_GENERATOR=STGZ -DBUILD_FORTRAN=ON .. \
    && make -j 4 \
    && make package \
    && chmod +x EnergyPlus-$ENERGYPLUS_VERSION-$ENERGYPLUS_SHA-Linux-x86_64.sh \
    && echo "y\r" | ./EnergyPlus-$ENERGYPLUS_VERSION-$ENERGYPLUS_SHA-Linux-x86_64.sh

# Remove the broken symlinks
#RUN cd /usr/local/bin \
#    && find -L . -type l -delete

# Add in the test files
#ADD test /usr/local/EnergyPlus-$ENERGYPLUS_INSTALL_VERSION/test_run
#RUN cp /usr/local/EnergyPlus-$ENERGYPLUS_INSTALL_VERSION/Energy+.idd \
#        /usr/local/EnergyPlus-$ENERGYPLUS_INSTALL_VERSION/test_run/

VOLUME /var/simdata
WORKDIR /var/simdata

CMD [ "/bin/bash" ]
