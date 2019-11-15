# -----------------------------------------------------------------------------
#
# challenge specific code and configuration
#
# -----------------------------------------------------------------------------
# Use rpy2 if you have R scoring functions
import rpy2.robjects as robjects
# import rpy2.robjects.packages as rpackages
# from rpy2.robjects.vectors import StrVector
import os
# # setup the R environment with tidyverse package:
# package_names = ("dplyr", "tidyr", "readr")
# if all(rpackages.isinstalled(x) for x in package_names):
#     have_package = True
# else:
#     have_package = False

# if not have_package:
#     utils = rpackages.importr('utils')
#     utils.chooseCRANmirror(ind=1)
#     packagenames_to_install = [x for x in package_names if not rpackages.isinstalled(x)]
#     if len(packagenames_to_install) > 0:
#         utils.install_packages(StrVector(packagenames_to_install))

# scoring scripts:
# __file__ = os.path.join(os.getcwd(), "challenge_config.py")
script_dir = os.path.dirname(os.path.abspath(__file__))
score1_path = os.path.join(script_dir, '../scoring_scripts', 'score_sc1.R')
score2_path = os.path.join(script_dir, '../scoring_scripts', 'score_sc2.R')
score3_path = os.path.join(script_dir, '../scoring_scripts', 'score_sc3.R')
score4_path = os.path.join(script_dir, '../scoring_scripts', 'score_sc4.R')

robjects.r("source('{}')".format(score1_path))
robjects.r("source('{}')".format(score2_path))
robjects.r("source('{}')".format(score3_path))
robjects.r("source('{}')".format(score4_path))

score_sc1 = robjects.r('score_sc1')
score_sc2 = robjects.r('score_sc2')
score_sc3 = robjects.r('score_sc3')
score_sc4 = robjects.r('score_sc4')

# validation scripts in R
validate1_path = os.path.join(
    script_dir, '../scoring_scripts', 'validate_sc1.R')
validate2_path = os.path.join(
    script_dir, '../scoring_scripts', 'validate_sc2.R')
validate3_path = os.path.join(
    script_dir, '../scoring_scripts', 'validate_sc3.R')
validate4_path = os.path.join(
    script_dir, '../scoring_scripts', 'validate_sc4.R')

robjects.r("source('{}')".format(validate1_path))
robjects.r("source('{}')".format(validate2_path))
robjects.r("source('{}')".format(validate3_path))
robjects.r("source('{}')".format(validate4_path))

validate_sc1 = robjects.r('validate_sc1')
validate_sc2 = robjects.r('validate_sc2')
validate_sc3 = robjects.r('validate_sc3')
validate_sc4 = robjects.r('validate_sc4')

# Configuring them here as a list will save a round-trip to the server
# every time the script starts and you can link the challenge queues to
# the correct scoring/validation functions.
# Predictions will be validated and scored


def validate_py_sc1(submission_path, goldstandard_path):
    '''
    Validate submission for subchallenge 1. 

    MUST USE ASSERTION ERRORS!!!

    eg.
    >>> assert os.path.basename(submission_path) == "prediction.tsv", \
    >>> "Submission file must be named prediction.tsv"
    or raise AssertionError()...
    Only assertion errors will be returned to participants,
    all other errors will be returned to the admin

    Args:
        submission_path:  Path to submission file
        goldstandard_path: Path to truth file

    Returns:
        Must return a boolean and validation message
    '''
    from synapseclient import Submission
    # Sometimes participants accidentally submit Projects/Folders
    assert not isinstance(submission_path, Submission), \
        "Submission must be a Synapse File and not Project/Folder"

    validation_result = validate_sc1(submission_path, goldstandard_path)
    assert validation_result[0][0] == 0, validation_result[1][0][:1000]

    is_valid = True
    message = "Passed Validation"
    return(is_valid, message)


def validate_py_sc2(submission_path, goldstandard_path):
    '''
    Validate submission for subchallenge 2. 

    MUST USE ASSERTION ERRORS!!!

    eg.
    >>> assert os.path.basename(submission_path) == "prediction.tsv", \
    >>> "Submission file must be named prediction.tsv"
    or raise AssertionError()...
    Only assertion errors will be returned to participants,
    all other errors will be returned to the admin

    Args:
        submission_path:  Path to submission file
        goldstandard_path: Path to truth file

    Returns:
        Must return a boolean and validation message
    '''
    from synapseclient import Submission
    # Sometimes participants accidentally submit Projects/Folders
    assert not isinstance(submission_path, Submission), \
        "Submission must be a Synapse File and not Project/Folder"

    validation_result = validate_sc2(submission_path, goldstandard_path)
    assert validation_result[0][0] == 0, validation_result[1][0][:1000]

    is_valid = True
    message = "Passed Validation"
    return(is_valid, message)


def validate_py_sc3(submission_path, goldstandard_path):
    '''
    Validate submission for subchallenge 3.

    MUST USE ASSERTION ERRORS!!!

    eg.
    >>> assert os.path.basename(submission_path) == "prediction.tsv", \
    >>> "Submission file must be named prediction.tsv"
    or raise AssertionError()...
    Only assertion errors will be returned to participants,
    all other errors will be returned to the admin

    Args:
        submission_path:  Path to submission file
        goldstandard_path: Path to truth file

    Returns:
        Must return a boolean and validation message
    '''
    from synapseclient import Submission
    # Sometimes participants accidentally submit Projects/Folders
    assert not isinstance(submission_path, Submission), \
        "Submission must be a Synapse File and not Project/Folder"

    validation_result = validate_sc3(submission_path, goldstandard_path)
    assert validation_result[0][0] == 0, validation_result[1][0][:1000]

    is_valid = True
    message = "Passed Validation"
    return(is_valid, message)


def validate_py_sc4(submission_path, goldstandard_path):
    '''
    Validate submission for subchallenge 4.

    MUST USE ASSERTION ERRORS!!!

    eg.
    >>> assert os.path.basename(submission_path) == "prediction.tsv", \
    >>> "Submission file must be named prediction.tsv"
    or raise AssertionError()...
    Only assertion errors will be returned to participants,
    all other errors will be returned to the admin

    Args:
        submission_path:  Path to submission file
        goldstandard_path: Path to truth file

    Returns:
        Must return a boolean and validation message
    '''
    from synapseclient import Submission
    # Sometimes participants accidentally submit Projects/Folders
    assert not isinstance(submission_path, Submission), \
        "Submission must be a Synapse File and not Project/Folder"

    validation_result = validate_sc4(submission_path, goldstandard_path)
    assert validation_result[0][0] == 0, validation_result[1][0][:1000]

    is_valid = True
    message = "Passed Validation"
    return(is_valid, message)


def validate_writeup(submission, goldstandard_path,
                     public=True, admin=None):
    '''
    Validates challenge writeup

    Args:
        submission: Submission object
        goldstandard_path: Unused
        syn: Synapse object
        public: If the writeup needs to be public. Defaults to True
        admin: Specify Synapse userid that writeup needs to be
               shared with
    Returns:
        (True, message) if validated, (False, message) if
        validation fails or throws exception
    '''
    from synapseclient.exceptions import SynapseHTTPError
    from synapseclient import AUTHENTICATED_USERS
    import synapseclient

    syn = synapseclient.login()
    # Add in users to share this with
    share_with = []
    try:
        if public:
            message = "Please make your private project ({}) public".format(
                submission['entityId'])
            share_with.append(message)
            ent = \
                syn.getPermissions(submission['entityId'], AUTHENTICATED_USERS)
            assert "READ" in ent and "DOWNLOAD" in ent, message
            ent = syn.getPermissions(submission['entityId'])
            assert "READ" in ent, message
        if admin is not None:
            message = (
                "Please share your private directory ({}) with the Synapse"
                " user `{}` with `Can Download` permissions.".format(
                    submission['entityId'], admin))
            share_with.append(message)
            ent = syn.getPermissions(submission['entityId'], admin)
            assert "READ" in ent and "DOWNLOAD" in ent, message
    except SynapseHTTPError as e:
        if e.response.status_code == 403:
            raise AssertionError("\n".join(share_with))
        else:
            raise(e)
    return True, "Validated!"


def score1(submission_path, goldstandard_path):
    '''
    Scoring function for subchallenge 1

    Args:
        submission_path:  Path to submission file
        goldstandard_path: Path to truth file

    Returns:
        Must return score dictionary and a scoring message
    '''
    sc = score_sc1(submission_path, goldstandard_path)
    score_dict = dict(score=round(sc[0], 6), round=3)
    score = score_dict.get('score')
    message = "Your submission has been scored!"
    return(score_dict, message)


def score2(submission_path, goldstandard_path):
    '''
    Scoring function for subchallenge 2

    Args:
        submission_path:  Path to submission file
        goldstandard_path: Path to truth file

    Returns:
        Must return score dictionary and a scoring message
    '''
    # Score against goldstandard
    sc = score_sc2(submission_path, goldstandard_path)

    score_dict = dict(score=round(sc[0], 6), round=3)
    score = score_dict.get('score')
    message = "Your submission has been scored!"
    return(score_dict, message)


def score3(submission_path, goldstandard_path):
    '''
    Scoring function for subchallenge 3

    Args:
        submission_path:  Path to submission file
        goldstandard_path: Path to truth file

    Returns:
        Must return score dictionary and a scoring message
    '''
    sc = score_sc3(submission_path, goldstandard_path)

    score_dict = dict(score=round(sc[0], 6), round=3)
    score = score_dict.get('score')
    message = "Your submission has been scored!"
    return(score_dict, message)


def score4(submission_path, goldstandard_path):
    '''
    Scoring function for subchallenge 4

    Args:
        submission_path:  Path to submission file
        goldstandard_path: Path to truth file

    Returns:
        Must return score dictionary and a scoring message
    '''
    # Score against goldstandard
    sc = score_sc4(submission_path, goldstandard_path)

    score_dict = dict(score=round(sc[0], 6), round=3)
    score = score_dict.get('score')
    message = "Your submission has been scored!"
    return(score_dict, message)


EVALUATION_QUEUES_CONFIG = [
    {
        'id': 9614297,
        'scoring_func': score1,
        'validation_func': validate_py_sc1,
        'goldstandard_path': os.path.join(script_dir, 'sc1gold.csv')
    },
    {
        'id': 9614298,
        'scoring_func': score2,
        'validation_func': validate_py_sc2,
        'goldstandard_path': os.path.join(script_dir, 'sc2gold.csv')
    },
    {
        'id': 9614302,
        'scoring_func': score3,
        'validation_func': validate_py_sc3,
        'goldstandard_path': os.path.join(script_dir, 'sc3gold.csv')
    },
    {
        'id': 9614303,
        'scoring_func': score4,
        'validation_func': validate_py_sc4,
        'goldstandard_path': os.path.join(script_dir, 'sc4gold.csv')
    },
    {
        'id': 9614264,
        'scoring_func': score1,  # will never run
        'validation_func': validate_writeup,
        'goldstandard_path': None
    }
]
