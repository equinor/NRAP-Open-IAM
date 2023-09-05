*******************
Introduction to git
*******************

.. toctree::

What is Git?
============

Git is one of the most often used version control systems in the world.
The main benefits and features provided for developers employing the version
control system are the following:

* system provides means for developers to keep track of code changes,

* it allows developers to see a history of changes, to work on the same code
  pieces at the same time, to isolate their code through branching,
  merge code from different branches,

* helps developers see and resolve conflicts on code merges, to revert
  their changes to a previous state, to merge only selected changes.

Git is an example of what is called a Distributed Version Control System:
each developer has a copy of the whole repository on their computer and
can see the entire history of changes. Git is also cross-platform which makes
it useful when developers are working with different OSs.
Developers can divide their work between different branches depending
on different priorities. Developers can create experimental features and
revert the changes when something goes wrong.

The NRAP-Open-IAM project is hosted at a public repository located at
https://gitlab.com/NRAP/OpenIAM. `GitLab <https://gitlab.com>`_ is a web-based
Git-repository manager with wiki and issue-tracking
features, using an open-source license. There is also the NRAP-Open-IAM development
repository which is private and open only to the members of the development group.
The master branch of the development repository NRAP-Open-IAM code
is protected. Any new code development is done on separate branches and goes through a
review process before merging with the master branch. It ensures that new code
does not break existing features and is properly documented and tested.
A merge request is approved when the code is complete and ready to be part
of the master branch.

In order to start working with the NRAP-Open-IAM the developer should register
on https://gitlab.com and provide the user name to be added. Once developer
has registered, the code can be downloaded in a compressed format from GitLab.
However, the compressed file will only contain the source code; the git repository
files including revision history will not be included. To utilize Git during
development, the repository should be cloned with the command ::

    git clone https://gitlab.com/NRAP/OpenIAM

Note that the command uses address for NRAP-Open-IAM repository open for public.
The command creates a directory named *OpenIAM* (at the current
local file system location), initializes a .git directory inside it,
pulls down all the data for that repository, and checks out a working
copy of the latest version.

If the clone of the repository is needed in a different location, e.g., folder
with name *new_NRAPOpenIAM_location*, then the following command should be used: ::

    git clone https://gitlab.com/NRAP/OpenIAM new_NRAPOpenIAM_location

This command does the same thing as the previous one, but the target directory
is now called *new_NRAPOpenIAM_location*.

The installation of Git on the local computer depends on the operating system.
The current NRAP-Open-IAM developers use the following Git management systems:

* Windows: Git installed with Cygwin, Git Bash

* Linux, Mac: Git might already be installed locally.

If Git is not available by default, it is worth checking
the following sources: for Mac https://git-scm.com/download/gui/mac, and
for Linux https://git-scm.com/download/gui/linux.
Additional Git clients compatible with Windows operating systems can be found at
https://git-scm.com/download/gui/windows .

Comprehensive book on Git which includes more details that can be covered
in this guide is available as an open-source project at https://git-scm.com/book/en/v2.
The website https://gitlab.com provides an introduction material on use of
GitLab and Git in general.


Most frequently used Git commands
=================================

As the name of the current section suggests we provide a list of the most
frequently used commands in Table 1.


.. table:: Most common Git commands
    :widths: 55 45

    +---------------------------------------------------+----------------------------------------------------+
    |                  **Command**                      |                         **Description**            |
    +===================================================+====================================================+
    | **Get Help**                                                                                           |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git help``                                    | Get help on git.                                   |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git help <command>``                          | Get help on any git command.                       |
    +                                                   |                                                    +
    | or                                                |                                                    |
    +                                                   |                                                    +
    | $ ``git <command> -h``                            |                                                    |
    +---------------------------------------------------+----------------------------------------------------+
    | **Configure**                                                                                          |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git config --global user.name <name>``        | Set the name the developer wants to be attached    |
    +                                                   |                                                    +
    |                                                   | to the commits.                                    |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git config --global user.email <email>``      | Set the email the developer wants to be attached   |
    +                                                   |                                                    +
    |                                                   | to the commits.                                    |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git config --list``                           | List git configuration settings.                   |
    +---------------------------------------------------+----------------------------------------------------+
    | **Create and Clone**                                                                                   |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git init <new-repository>``                   | Create a new local repository.                     |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git clone <url>``                             | Clone an existing repository.                      |
    +---------------------------------------------------+----------------------------------------------------+
    | **Keep Track of Changes**                                                                              |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git status``                                  | List all new and modified files in the working     |
    +                                                   |                                                    +
    |                                                   | directory.                                         |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git diff``                                    | Show not yet staged differences in the tracked     |
    +                                                   |                                                    +
    |                                                   | files.                                             |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git diff --staged``                           | Show differences between staged and committed      |
    +                                                   |                                                    +
    |                                                   | version of files.                                  |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git add <file or dir>``                       | Stage <file> or <dir> (take a snapshot             |
    +                                                   |                                                    +
    |                                                   | of the content) in preparation for commit          |
    +                                                   |                                                    +
    |                                                   | or stash.                                          |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git add -p <file>``                           | Stage selected changes in <file>.                  |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git add -A``                                  | Stage all changes (new untracked, modified         |
    +                                                   |                                                    +
    |                                                   | tracked and deleted files) in the working          |
    +                                                   |                                                    +
    |                                                   | directory.                                         |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git add .``                                   | Stage all new untracked and modified tracked       |
    +                                                   |                                                    +
    |                                                   | (but not deleted files) in the working directory.  |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git add -u``                                  | Stage all modified tracked and not deleted files   |
    +                                                   |                                                    +
    |                                                   | in the working directory.                          |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git reset HEAD <file>``                       | Unstage <file> but preserve its content.           |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git checkout -- <file>``                      | Revert the convert of <file> to the last commit    |
    +                                                   |                                                    +
    |                                                   | version.                                           |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git commit -m <commit-message>``              | Record (commit) previously staged changes in       |
    +                                                   |                                                    +
    |                                                   | version history.                                   |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git commit -a``                               | Commit all changes in tracked files (without       |
    +                                                   |                                                    +
    |                                                   | staging step).                                     |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git commit --amend``                          | Change the last commit.                            |
    +---------------------------------------------------+----------------------------------------------------+
    | **Review Commit History**                                                                              |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git log``                                     | Show all commits starting with the most recent.    |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git log -p <file>``                           | Show changes over time for a given <file>.         |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git blame <file>``                            | Show who changed what and when in a given <file>.  |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git diff <branch1>...<branch2>``              | Show differences in the content of two branches.   |
    +---------------------------------------------------+----------------------------------------------------+
    | **Create a Branch**                                                                                    |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git branch``                                  | Show all local branches in the current repository. |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git branch -av``                              | Show all existing branches with some extra         |
    +                                                   |                                                    +
    |                                                   | information.                                       |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git branch <new-branch-name>``                | Create a new branch based on the current branch.   |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git branch -d <branch>``                      | Delete fully merged <branch>.                      |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git checkout <branch>``                       | Switch to different <branch>.                      |
    +---------------------------------------------------+----------------------------------------------------+
    | **Update**                                                                                             |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git remote -v``                               | Show all currently configured remotes.             |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git remote show <name-of-remote>``            | Show information about a remote.                   |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git fetch <name-of-remote>``                  | Download all changes from the                      |
    +                                                   |                                                    +
    |                                                   | <name-of-remote> branch but do not integrate       |
    +                                                   |                                                    +
    |                                                   | them into the local branch.                        |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git pull <name-of-remote> <branch>``          | Download all changes from the                      |
    +                                                   |                                                    +
    |                                                   | <name-of-remote> branch and integrate/merge        |
    +                                                   |                                                    +
    |                                                   | them into the local <branch>.                      |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git push <name-of-remote> <branch>``          | Publish local changes from <branch> to the         |
    +                                                   |                                                    +
    |                                                   | <name-of-remote>.                                  |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git merge <branch>``                          | Merge <branch> to the current local branch.        |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git rebase <branch>``                         | Rebase the current local branch onto <branch>.     |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git rebase --continue``                       | Continue a rebase after resolving conflicts.       |
    +---------------------------------------------------+----------------------------------------------------+
    | **Save Fragments**                                                                                     |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git stash``                                   | Record the current state of the working directory  |
    +                                                   |                                                    +
    | or                                                | and go back to the clean working directory.        |
    +                                                   |                                                    +
    | $ ``git stash -m <stash-entry-message>``          |                                                    |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git stash apply``                             | Restore the most recently stashed files.           |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git stash list``                              | List the currently available stash entries.        |
    +---------------------------------------------------+----------------------------------------------------+
    | $ ``git stash clear``                             | Remove all the stash entries.                      |
    +---------------------------------------------------+----------------------------------------------------+
