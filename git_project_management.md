# 2026-ysjin-seqdep

## This branch is hive!!

### common git commands
```bash
git branch --show current # check which branch you are --show-current

git remote -v # check the remote repo urls for fetch and pull
git remote set-url --push origin DISABLE # if you want to disable git push for a branch
git remote set-url --push origin <real-git-repo-url> # if you want to reset git push for a branch

# The following commands will contaminate main branch
# 1) switch to main and push
git checkout main
git push
# 2) merge hive branch to main
git checkout main
git merge hive
git push
# 3) force push and overwrite main branch
git push origin hive:main

# Dry run command
git push --dry-run # Would push to origin hive-real-data, confirm before running git push

# Tags
# when you edited codes, if this is important for activities like paper revision & plotting, use tags to track the codes/commits.
# example usage:
git status
git add .
git commit -m "added R script for plotting"
git tag -a run_20260302 -m "generate figure 1"
git push --tags
# if you want to check tags
git tag
git show run_20260302 # show the description of a tag
git checkout run_20260302 # switch to the version for run_20260302 tag

# dev on main
git checkout main
git pull origin main

git add -A
git commit -m "method update"
git push

# run on hive
git checkout hive
git fetch origin
git merge origin/main
git add -A
git commit -m "hive: small adjustments for real run"
git push origin hive # if this is the first push, git push --set-upstream origin hive
```
