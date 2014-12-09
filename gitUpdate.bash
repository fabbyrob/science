#usage:
echo 'bash gitUpdate.bash "/Document/repos/git/science/" myChangedFile.txt "my awesome comment about my change"'

cd $1 #change to git repo
pwd

echo "pulling..."
git pull 
echo "adding..."
git add $2 #add the file
echo "committing..."
while [ -n "$3" ]; do
	message = "$message $3"
	shift
done
git commit -m \'$message\' #commit changes, with a comment
git push origin master
