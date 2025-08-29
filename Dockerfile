# Base image
FROM rocker/shiny:4.4.1

# General updates
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git libxml2-dev libmagick++-dev libssl-dev libharfbuzz-dev libfribidi-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# --- Add these lines: ensure C++14 is set for package builds ---
RUN echo "CXX14FLAGS=-std=c++14" >> /usr/local/lib/R/etc/Makeconf
RUN echo "CXX14 = g++" >> /usr/local/lib/R/etc/Makeconf
RUN echo "CXX14STD = -std=c++14" >> /usr/local/lib/R/etc/Makeconf
# --------------------------------------------------------------

# Install the required packages
RUN Rscript -e 'install.packages(c("renv"))'
COPY /renv.lock /srv/shiny-server/renv.lock
RUN Rscript -e 'setwd("/srv/shiny-server/");renv::restore();'

# Copy the app files (scripts, data, etc.)
RUN rm -rf /srv/shiny-server/*
COPY /app/ /srv/shiny-server/

# Ensure that the expected user is present in the container
RUN if id shiny &>/dev/null && [ "$(id -u shiny)" -ne 999 ]; then \
        userdel -r shiny; \
        id -u 999 &>/dev/null && userdel -r "$(id -un 999)"; \
    fi; \
    useradd -u 999 -m -s /bin/bash shiny; \
    chown -R shiny:shiny /srv/shiny-server/ /var/lib/shiny-server/ /var/log/shiny-server/

USER shiny
EXPOSE 3838
CMD ["/usr/bin/shiny-server"]
