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
# Authors:  Brian D. Weitzner (brian.weitzner@gmail.com)                      #
#           Tim Jacobs (TimJacobs2@gmail.com)                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if [ -z "$1" ]; then
	echo "Your GitHub user name must be the first argument to this script!"
	exit
fi

if [ -z "$2" ]; then
    echo "Your GitHub password must be the second argument to this script!"
fi

github_user_name=$1
github_password=$2
repo="rosetta"

read -p "Where would you like to clone $repo? " path
if [ -z $path ]; then
	path="."
fi

if [ ! -d $path ]; then
	echo "\033[0;33m'$path' does not exist!\033[0m You'll need to create '$path' if you want to install rosetta there."
	while true; do
		read -p "Would you like to create this directory now? " yn
		case $yn in
			[Yy] | [Yy][Ee][Ss] ) mkdir $path; break;;
			[Nn] | [Nn][Oo] ) exit;;
		* ) echo "Please answer yes (y) or no (n).";;
		esac
	done
fi  

while true; do
	read -p "Would you like to clone over ssh (s) or https (h) - Note that ssh keys are required for cloning over ssh (Default: ssh)? " protocol
	case $protocol in
		[Ss] | [Ss][Ss][Hh] | "" ) url=git@github.com:RosettaCommons/; break;;
		[Hh] | [Hh][Tt][Tt][Pp][Ss] ) url=https://$github_user_name:$github_password@github.com/RosettaCommons/; break;;
	*) echo "Please answer ssh (s) or https (h).";;
	esac
done

path="$path/"

echo "\033[0;34mCloning Rosetta...\033[0m"
hash git >/dev/null && /usr/bin/env git clone $url$repo.git $path$repo || {
echo "Can't clone! It's likely that git is not installed and/or you are cloning over SSH without ssh keys setup.\nSee https://help.github.com/articles/error-permission-denied-publickey for instructions on how to setup SSH keys for github."
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

starting_dir=$PWD
cd $path/$repo

echo "\033[0;34mDisabling fast-forward merges on master...\033[0m"
git config branch.master.mergeoptions "--no-ff"

echo "\033[0;34mConfiguring commit message template...\033[0m"
git config commit.template .commit_template.txt

cd .git/hooks
url="https://github.com/RosettaCommons/rosetta_tools/raw/master"
for hook in pre-commit post-commit; do 
	echo "\033[0;34mConfiguring the $hook hook...\033[0m"
	curl -u $github_user_name:$github_password -L $url/git_hooks/$hook > $hook
	chmod +x $hook
done

echo "\033[0;34mConfiguring aliases...\033[0m"
git config alias.tracked-branch '!sh -c "git checkout -b $1 && git push origin $1:$2/$1 && git branch --set-upstream $1  origin/$2/$1" -'
git config alias.personal-tracked-branch '!sh -c "git tracked-branch $1 $github_user_name" -'
sed -ie "s/\$github_user_name/$github_user_name/g" .git/config

git config --global alias.show-graph "log --graph --abbrev-commit --pretty=oneline"

echo "\033[0;34mConfiguring git colors...\033[0m"
git config --global color.branch auto
git config --global color.diff auto
git config --global color.interactive auto
git config --global color.status auto

cd $starting_dir
echo "\033[0;32mDone configuring your Rosetta git repository!\033[0m"
echo "\033[0;32mRemember to check out rosetta_tools and rosetta_demos repositories!\033[0m"
