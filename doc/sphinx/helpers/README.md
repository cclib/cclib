# Documentation helpers for developers

This directory contains helpers for working with documentation.

## [milestone_issues_and_prs.py](milestone_issues_and_prs.py)

The Sphinx-based changelog consists of the reduction of issues and PRs to the repository between each version tag.
This script generates a Sphinx fragment containing the raw issue and PR titles, numbers, and authors that have been assigned to a GitHub milestone.
These should then be summarized and combined into more useful changelog entries.

In order to run the script (`python milestone_issues_and_prs.py > milestone_issues_and_prs.rst`),

- export the environment variable `GITHUB_API_TOKEN` to contain a Personal Access Token that has ~read privileges
- change the name of the milestone in [milestone_issues_and_prs.graphql](milestone_issues_and_prs.graphql) to match the upcoming release milestone
