# keyfiles
Important Notes on Skills Learned

Git skills   
On resolving conflicts between branches:
I worked on branch "mengp2" forked from "master" branch in spring 2016. "mengp2" branch stopped updating since May 20ish 2016. 
Over the summer and fall, Dener started "dev" branch from the "master" branch and has accumulated lots of updates till now.
Now I want to resume my work on "mengp2" branch, but base it on top of the latest dev head.

There's not even dev branch exists in my local repository
git fetch               % Fetch branches from one or more other repositories; dev branch shows up
git checkout mengp2
git rebase dev          % cannot proceed due to conflicts of changes to the same section in the same file from the two branches
git rebase -Xours dev   % in favor of the dev branch,  X codes on mengp2 branch

Now my local mengp2 branch sits on top of latest dev branch, but is different from the remote mengp2 branch
git push origin mengp2        % cannot proceed due to conflicts between local and remote mengp2 branch
git push -f origin mengp2     % force the local revision to the remote repo 

Now problem solved! 
