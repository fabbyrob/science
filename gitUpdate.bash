#usage:
echo 'bash gitUpdate.bash "/Document/repos/git/science/" myChangedFile.txt "my awesome comment about my change"'

cd $1 #change to git repo
pwd

#git pull 
git add $1$2 #add the file
#git commit -m $3 #commit changes, with a comment
#git push origin master
