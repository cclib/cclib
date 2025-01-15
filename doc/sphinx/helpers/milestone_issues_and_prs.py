# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
import json
import sys
from pathlib import Path

from github_graphql_query import FILTER_AUTHORS, execute_query, transform_author

if __name__ == "__main__":
    raw = execute_query(Path("milestone_issues_and_prs.graphql"))
    result = json.loads(raw)

    # useful for debugging
    Path("result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    milestones = result["data"]["repository"]["milestones"]["nodes"]
    assert len(milestones) == 1
    milestone = milestones[0]

    issues = milestone["issues"]["nodes"]
    pull_requests = milestone["pullRequests"]["nodes"]

    # flatten the author field
    for issue in issues:
        transform_author(issue)
    for pull_request in pull_requests:
        transform_author(pull_request)
        for closing_issue in pull_request["closingIssuesReferences"]["nodes"]:
            transform_author(closing_issue)
        # collect all authors from commits
        all_authors = set()
        commits = [item["commit"] for item in pull_request["commits"]["nodes"]]
        for commit in commits:
            authors = commit["authors"]["nodes"]
            for author in authors:
                all_authors.add(author["user"]["login"])
        all_authors.add(pull_request["author"])
        pull_request["authors"] = all_authors - FILTER_AUTHORS
        # flatten nodes
        pull_request["closingIssuesReferences"] = pull_request["closingIssuesReferences"]["nodes"]

    all_closing_issues = set()

    for pull_request in pull_requests:
        for closing_issue in pull_request["closingIssuesReferences"]:
            if closing_issue not in issues:
                print(
                    f"PR {pull_request['url']} contains something not in issue list: {closing_issue}",
                    file=sys.stderr,
                )
            # all_closing_issues.add(tuple(sorted(closing_issue.items())))

    all_issues = set(tuple(tuple(sorted(x.items())) for x in issues))
    issues_without_pr = all_issues - all_closing_issues

    lines = list()
    lines.append("Issues")
    for issue in issues:
        lines.append(f"    * {issue['title']} (#{issue['number']})")
    lines.append("Pull requests")
    for pull_request in pull_requests:
        if pull_request["author"] not in FILTER_AUTHORS:
            authors = pull_request["authors"]
            lines.append(
                f"    * {pull_request['title']} ({', '.join(sorted(authors))}, #{pull_request['number']})"
            )
            closing_issues = pull_request["closingIssuesReferences"]
            for issue in closing_issues:
                lines.append(f"        * {issue['title']} (#{issue['number']})")

    print("\n".join(line.replace("`", "``") for line in lines))
