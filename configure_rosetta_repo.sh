# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available
# (c) under license.
# (c) The Rosetta software is developed by the contributing members of the
# (c) Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington UW
# (c) TechTransfer, email: license@u.washington.edu.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Brief:   This shell script clones repositories from GitHub and configures   #
#          them to play nicely with how our community is organized.           #
#                                                                             #
# Author:  Brian D. Weitzner (brian.weitzner@gmail.com)                       #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

repo="rosetta"

while true; do
    read -p "Where do you to clone $repo? " path
    if [ ! -d $path ]
    then
        echo "\033[0;33m\'$path\' does not exist!\033[0m You'll need to create $path if you want to install rosetta there."
        while true; do
            read -p "\033[0;34mWould you like to create this directory now?\033[0m " yn
            case $yn in
                [Yy]* ) mkdir $path; break;;
                [Nn]* ) exit;;
                * ) echo "Please answer yes or no.";;
            esac
        done
    fi  
    
    echo "\033[0;34mCloning Rosetta...\033[0m"
    hash git >/dev/null && /usr/bin/env git clone git@github.com:RosettaCommons/$repo.git $path/$repo || {
    echo "git is not installed!"
    exit
    }
                                    
    echo "\033[0;32m"'     ___           ___           ___           ___           __         ___         ___      '"\033[0m"
    echo "\033[0;32m"'    /\  \         /\  \         /\__\         /\__\         /\__\      /\__\       /\  \     '"\033[0m"
    echo "\033[0;32m"'   /::\  \       /::\  \       /:/ _/_       /:/ _/_       /:/  /     /:/  /      /::\  \    '"\033[0m"
    echo "\033[0;32m"'  /:/\:\__\     /:/\:\  \     /:/ /\  \     /:/ /\__\     /:/__/     /:/__/      /:/\:\  \   '"\033[0m"
    echo "\033[0;32m"' /:/ /:/  /    /:/  \:\  \   /:/ /::\  \   /:/ /:/ _/_   /::\  \    /::\  \     /:/ /::\  \  '"\033[0m"
    echo "\033[0;32m"'/:/_/:/__/___ /:/__/ \:\__\ /:/_/:/\:\__\ /:/_/:/ /\__\ /:/\:\  \  /:/\:\  \   /:/_/:/\:\__\ '"\033[0m"
    echo "\033[0;32m"'\:\/:::::/  / \:\  \ /:/  / \:\/:/ /:/  / \:\/:/ /:/  / \/__\:\  \ \/__\:\  \  \:\/:/  \/__/ '"\033[0m"
    echo "\033[0;32m"' \::/~~/~~~~   \:\  /:/  /   \::/ /:/  /   \::/_/:/  /       \:\  \     \:\  \  \::/__/      '"\033[0m"
    echo "\033[0;32m"'  \:\~~\        \:\/:/  /     \/_/:/  /     \:\/:/  /         \:\  \     \:\  \  \:\  \      '"\033[0m"
    echo "\033[0;32m"'   \:\__\        \::/  /        /:/  /       \::/  /           \:\__\     \:\__\  \:\__\     '"\033[0m"
    echo "\033[0;32m"'    \/__/         \/__/         \/__/         \/__/             \/__/      \/__/   \/__/     '"\033[0m"

    echo "\n\n \033[0;32m....is now cloned.\033[0m"

    break;
done
