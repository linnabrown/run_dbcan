"""A Sphinx extension for linking to your project's issue tracker."""
import re
from typing import Callable, Optional, Tuple

from docutils import nodes, utils
from sphinx.config import Config
from sphinx.util.nodes import split_explicit_title

__version__ = "3.0.1"
__author__ = "Steven Loria"
__license__ = "MIT"


def cve_role(name, rawtext, text, lineno, inliner, options=None, content=None):
    """Sphinx role for linking to a CVE on https://cve.mitre.org.

    Examples: ::

        :cve:`CVE-2018-17175`

    """
    options = options or {}
    content = content or []
    has_explicit_title, title, target = split_explicit_title(text)

    target = utils.unescape(target).strip()
    title = utils.unescape(title).strip()
    ref = f"https://cve.mitre.org/cgi-bin/cvename.cgi?name={target}"
    text = title if has_explicit_title else target
    link = nodes.reference(text=text, refuri=ref, **options)
    return [link], []


GITHUB_USER_RE = re.compile("^https://github.com/([^/]+)/([^/]+)/.*")


def _get_default_group_and_project(config: Config, uri_config_option: str) -> Optional[Tuple[str, str]]:
    """
    Retrieves the default group and project names from the configuration.

    This function extracts and returns the default group and project names based on the configuration settings. It supports both legacy and new configuration options, raising an error if both are defined.

    Parameters
    ----------
        config (Config): The configuration object containing the settings.
        uri_config_option (str): The URI configuration option name.

    Returns
    -------
        Optional[Tuple[str, str]]: A tuple containing the group and project names, or None if not set.

    Raises
    ------
        ValueError: If both old and new configuration options are set, or if the group/project format is incorrect.
    """
    old_config = getattr(config, "issues_github_path", None)
    new_config = getattr(config, "issues_default_group_project", None)

    if old_config and new_config:
        raise ValueError(
            "Both 'issues_github_path' and 'issues_default_group_project' are set, even"
            " though they define the same setting.  "
            "Please only define one of these."
        )
    group_and_project = new_config or old_config

    if group_and_project:
        assert isinstance(group_and_project, str)
        try:
            group, project = group_and_project.split("/", maxsplit=1)
            return group, project
        except ValueError as e:
            raise ValueError(
                "`issues_github_path` or `issues_default_group_project` needs to "
                "define a value in the form of `<group or user>/<project>` "
                f"but `{config}` was given."
            ) from e

    # If group and project was not set, we need to look for it within the github url
    # for backward compatibility
    if not group_and_project:
        uri = getattr(config, uri_config_option)
        if uri:
            match = GITHUB_USER_RE.match(uri)
            if match:
                return match.groups()[0], match.groups()[1]

    return None


def _get_placeholder(uri_config_option: str) -> str:
    """
    Extracts a placeholder value from a URI configuration option.

    This function processes a URI configuration option name to extract a meaningful placeholder. The placeholder is typically a part of the configuration option name.

    Parameters
    ----------
        uri_config_option (str): The URI configuration option name.

    Returns
    -------
        str: The extracted placeholder string.

    Note:
        The function handles different naming conventions for the configuration options.
    """
    try:
        # i.e. issues_pr_uri -> pr
        return uri_config_option[:-4].split("_", maxsplit=1)[1]
    except IndexError:
        # issues_uri -> issue
        return uri_config_option[:-5]


def _get_uri_template(
    config: Config,
    uri_config_option: str,
) -> str:
    """
    Get a URL format template that can be filled with user information based on the given configuration

    The result always contains the following placeholder
      - n (the issue number, user, pull request, etc...)

    The result can contain the following other placeholders
      - group (same as user in github)
      - project

    Examples for possible results:

         - "https://github.com/{group}/{project}/issues/{n}"

         - "https://gitlab.company.com/{group}/{project}/{n}"

         - "https://fancy.issuetrack.com?group={group}&project={project}&issue={n}"

    Raises
    ------
         - ValueError if the given uri contains an invalid placeholder
    """
    format_string = str(getattr(config, uri_config_option))
    placeholder = _get_placeholder(uri_config_option)

    result = format_string.replace(f"{{{placeholder}}}", "{n}")

    try:
        result.format(project="", group="", n="")
    except (NameError, KeyError) as e:
        raise ValueError(
            f"The `{uri_config_option}` option contains invalid placeholders. "
            f"Only {{group}}, {{projects}} and {{{placeholder}}} are allowed."
            f'Invalid format string: "{format_string}".'
        ) from e
    return result


def _get_uri(
    uri_config_option: str,
    config: Config,
    number: str,
    group_and_project: Optional[Tuple[str, str]] = None,
) -> str:
    """
    Constructs a URI based on configuration options and provided parameters.

    This function generates a URI using a format string obtained from the configuration, replacing placeholders with actual values provided as parameters or from the configuration. It supports backward compatibility by allowing replacement of default group/project in the format string.

    Parameters
    ----------
        uri_config_option (str): The configuration option name that specifies the URI format.
        config (Config): The configuration object containing settings and format strings.
        number (str): A string, typically a number, to be included in the URI.
        group_and_project (Optional[Tuple[str, str]]): A tuple containing the group and project names. If not provided, the default from the configuration is used.

    Returns
    -------
        str: The constructed URI based on the provided information and configuration.

    Raises
    ------
        ValueError: If the format string requires a group/project to be defined, and it is not provided in the function call or configuration.
    """
    format_string = _get_uri_template(config, uri_config_option)

    url_vars = {"n": number}

    config_group_and_project = _get_default_group_and_project(config, uri_config_option)
    if group_and_project:
        # Group and Project defined by call
        if config_group_and_project:
            to_replace = "/".join(config_group_and_project)
            if to_replace in format_string:
                # Backward compatibility, replace default group/project
                # with {group}/{project}
                format_string = format_string.replace(to_replace, "{group}/{project}")
        (url_vars["group"], url_vars["project"]) = group_and_project
    elif config_group_and_project:
        # If not defined by call use the default if given
        (url_vars["group"], url_vars["project"]) = config_group_and_project

    try:
        return format_string.format(**url_vars)
    except (NameError, KeyError) as e:
        # The format string was checked before, that it contains no additional not
        # supported placeholders. So this occur
        raise ValueError(
            f"The `{uri_config_option}` format `{format_string}` requires a "
            f"group/project to be defined in `issues_default_group_project`."
        ) from e


def cwe_role(name, rawtext, text, lineno, inliner, options=None, content=None):
    """Sphinx role for linking to a CWE on https://cwe.mitre.org.

    Examples: ::

        :cwe:`CWE-787`

    """
    options = options or {}
    content = content or []
    has_explicit_title, title, target = split_explicit_title(text)

    target = utils.unescape(target).strip()
    title = utils.unescape(title).strip()
    number = target[4:]
    ref = f"https://cwe.mitre.org/data/definitions/{number}.html"
    text = title if has_explicit_title else target
    link = nodes.reference(text=text, refuri=ref, **options)
    return [link], []


class IssueRole:
    """
    A class for formatting and linking issues, pull requests, merge requests, and commits.

    This class handles the generation of links to issues, pull requests, merge requests, and commits based on a configuration prefix and optionally provided text formatting. It supports both internal and external repository references.

    Attributes
    ----------
        ELEMENT_SEPARATORS (str): Symbols used to separate elements in references.
        EXTERNAL_REPO_REGEX (re.Pattern): Regular expression for matching external repository references.

    Methods
    -------
        __init__: Initializes the IssueRole instance.
        default_pre_format_text: Default text formatting method.
        format_text: Formats text with supported separators.
        make_node: Creates a docutils node for the reference.
        __call__: Processes text and returns formatted reference nodes.
    """

    # Symbols used to separate and issue/pull request/merge request etc
    # i.e
    #   - group/project#2323 for issues
    #   - group/project!1234 for merge requests (in gitlab)
    #   - group/project@adbc1234 for commits
    ELEMENT_SEPARATORS = "#@!"

    EXTERNAL_REPO_REGEX = re.compile(rf"^(\w+)/(.+)([{ELEMENT_SEPARATORS}])([\w]+)$")

    def __init__(
        self,
        config_prefix: str,
        pre_format_text: Callable[[Config, str], str] = None,
    ):
        """
        Initializes the IssueRole instance.

        Parameters
        ----------
            config_prefix (str): The prefix used for configuration options.
            pre_format_text (Callable[[Config, str], str], optional): A function for pre-formatting text before generating the reference.
        """
        self.uri_config = f"{config_prefix}_uri"
        self.separator_config = f"{config_prefix}_prefix"
        self.pre_format_text = pre_format_text or self.default_pre_format_text

    @staticmethod
    def default_pre_format_text(config: Config, text: str) -> str:
        """
        Default method for pre-formatting text.

        Parameters
        ----------
            config (Config): The configuration object.
            text (str): The text to format.

        Returns
        -------
            str: The formatted text.
        """
        return text

    def format_text(self, config: Config, issue_no: str) -> str:
        """
        Formats the issue number with the appropriate separator.

        Parameters
        ----------
            config (Config): The configuration object.
            issue_no (str): The issue number to format.

        Returns
        -------
            str: The formatted issue number with the separator.

        Raises
        ------
            ValueError: If an invalid separator is specified in the configuration.
        """
        separator = getattr(config, self.separator_config)
        if separator not in self.ELEMENT_SEPARATORS:
            raise ValueError(
                f"Option {self.separator_config} has to be one of " f"{', '.join(self.ELEMENT_SEPARATORS)}."
            )
        text = self.pre_format_text(config, issue_no.lstrip(self.ELEMENT_SEPARATORS))
        return f"{separator}{text}"

    def make_node(self, name: str, issue_no: str, config: Config, options=None):
        """
        Creates a docutils node for the given issue number.

        Parameters
        ----------
            name (str): The role name.
            issue_no (str): The issue number.
            config (Config): The configuration object.
            options (dict, optional): Additional options for the node.

        Returns
        -------
            docutils.nodes.reference: A reference node pointing to the issue.

        Note:
            Handles both internal and external repository references.
        """
        if issue_no in ("-", "0"):
            return None

        options = options or {}

        has_explicit_title, title, target = split_explicit_title(issue_no)

        if has_explicit_title:
            issue_no = str(target)

        repo_match = self.EXTERNAL_REPO_REGEX.match(issue_no)

        if repo_match:
            # External repo
            group, project, original_separator, issue_no = repo_match.groups()
            text = f"{group}/{project}{self.format_text(config, issue_no)}"
            ref = _get_uri(
                self.uri_config,
                config,
                issue_no,
                (group, project),
            )
        else:
            text = self.format_text(config, issue_no)
            ref = _get_uri(self.uri_config, config, issue_no)
        if has_explicit_title:
            return nodes.reference(text=title, refuri=ref, **options)
        else:
            return nodes.reference(text=text, refuri=ref, **options)

    def __call__(self, name, rawtext, text, lineno, inliner, options=None, content=None):
        """
        Processes the raw text and returns a list of reference nodes.

        Called by docutils when the role is invoked in the documentation.

        Parameters
        ----------
            name (str): The role name.
            rawtext (str): The entire markup snippet, including the role.
            text (str): The text marked with the role.
            lineno (int): The line number where the role occurs.
            inliner (docutils.parsers.rst.states.Inliner): The inliner instance.
            options (dict, optional): Directive options for further customization.
            content (list, optional): The directive content for nested parsing.

        Returns
        -------
            tuple: A two-item tuple containing a list of nodes and a list of system messages.
        """
        options = options or {}
        content = content or []
        issue_nos = [each.strip() for each in utils.unescape(text).split(",")]
        config = inliner.document.settings.env.app.config
        ret = []
        for i, issue_no in enumerate(issue_nos):
            node = self.make_node(name, issue_no, config, options=options)
            ret.append(node)
            if i != len(issue_nos) - 1:
                sep = nodes.raw(text=", ", format="html")
                ret.append(sep)
        return ret, []


"""Sphinx role for linking to an issue. Must have
`issues_uri` or `issues_default_group_project` configured in ``conf.py``.
Examples: ::
    :issue:`123`
    :issue:`42,45`
    :issue:`sloria/konch#123`
"""
issue_role = IssueRole(
    config_prefix="issues",
)

"""Sphinx role for linking to a pull request. Must have
`issues_pr_uri` or `issues_default_group_project` configured in ``conf.py``.
Examples: ::
    :pr:`123`
    :pr:`42,45`
    :pr:`sloria/konch#43`
"""
pr_role = IssueRole(
    config_prefix="issues_pr",
)


def format_commit_text(config, sha):
    """
    Formats a commit SHA to a shorter version.

    This function truncates the given commit SHA to its first 7 characters, which is a common short form for representing commit hashes.

    Parameters
    ----------
        config (Config): The configuration object. Currently not used in the function but included for consistency and potential future use.
        sha (str): The full commit SHA string.

    Returns
    -------
        str: A truncated version of the commit SHA, consisting of the first 7 characters.
    """
    return sha[:7]


"""Sphinx role for linking to a commit. Must have
`issues_commit_uri` or `issues_default_group_project` configured in ``conf.py``.
Examples: ::
    :commit:`123abc456def`
    :commit:`sloria/konch@123abc456def`
"""
commit_role = IssueRole(
    config_prefix="issues_commit",
    pre_format_text=format_commit_text,
)

"""Sphinx role for linking to a user profile. Defaults to linking to
GitHub profiles, but the profile URIS can be configured via the
``issues_user_uri`` config value.

Examples: ::

    :user:`sloria`

Anchor text also works: ::

    :user:`Steven Loria <sloria>`
"""
user_role = IssueRole(config_prefix="issues_user")


def setup(app):
    """
    Configures the Sphinx application with custom settings for issue tracking and formatting.

    This function is used to set up various configurations for linking issues, pull requests (PRs), commits, user profiles, and more in Sphinx documentation. It defines custom URI templates and prefixes for these entities and registers several roles for inline markup in reStructuredText.

    Parameters
    ----------
        app (sphinx.application.Sphinx): The Sphinx application object.

    Returns
    -------
        dict: A dictionary containing the extension version and compatibility flags for parallel read and write operations.

    Note:
        This function adds several configuration values to the Sphinx app, related to issue tracking and referencing in documentation. It also registers custom roles like 'issue', 'pr', 'user', 'commit', etc., for inline linking in the documentation.
    """
    # Format template for issues URI
    # e.g. 'https://github.com/sloria/marshmallow/issues/{issue}
    app.add_config_value(
        "issues_uri",
        default="https://github.com/{group}/{project}/issues/{issue}",
        rebuild="html",
        types=[str],
    )
    app.add_config_value("issues_prefix", default="#", rebuild="html", types=[str])
    # Format template for PR URI
    # e.g. 'https://github.com/sloria/marshmallow/pull/{issue}
    app.add_config_value(
        "issues_pr_uri",
        default="https://github.com/{group}/{project}/pull/{pr}",
        rebuild="html",
        types=[str],
    )
    app.add_config_value("issues_pr_prefix", default="#", rebuild="html", types=[str])
    # Format template for commit URI
    # e.g. 'https://github.com/sloria/marshmallow/commits/{commit}
    app.add_config_value(
        "issues_commit_uri",
        default="https://github.com/{group}/{project}/commit/{commit}",
        rebuild="html",
        types=[str],
    )
    app.add_config_value("issues_commit_prefix", default="@", rebuild="html", types=[str])
    # There is no seperator config as a format_text function is given

    # Default User (Group)/Project eg. 'sloria/marshmallow'
    # Called github as the package was working with github only before
    app.add_config_value("issues_github_path", default=None, rebuild="html", types=[str])
    # Same as above but with new naming to reflect the new functionality
    # Only on of both can be set
    app.add_config_value("issues_default_group_project", default=None, rebuild="html", types=[str])
    # Format template for user profile URI
    # e.g. 'https://github.com/{user}'
    app.add_config_value(
        "issues_user_uri",
        default="https://github.com/{user}",
        rebuild="html",
        types=[str],
    )
    app.add_config_value("issues_user_prefix", default="@", rebuild="html", types=[str])
    app.add_role("issue", issue_role)
    app.add_role("pr", pr_role)
    app.add_role("user", user_role)
    app.add_role("commit", commit_role)
    app.add_role("cve", cve_role)
    app.add_role("cwe", cwe_role)
    return {
        "version": __version__,
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
