FROM ubuntu:14.04

MAINTAINER Jason W. DeGraw jason.degraw@nrel.gov

# Downloading from Github
ENV ENERGYPLUS_REPO_URL https://github.com/NREL/EnergyPlus.git

# Set up, clone the repo, and build
RUN apt-get update && apt-get install -y git cmake g++ python \
    && mkdir energyplus \
    && cd energyplus \
    && git clone -b develop --single-branch $ENERGYPLUS_REPO_URL . \
    && mkdir build \
    && cd build \
    && cmake -DBUILD_TESTING=ON .. \
    && make -j 4 \
    && make install

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
