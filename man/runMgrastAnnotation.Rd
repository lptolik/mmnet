\name{runMgrastAnnotation}
\alias{runMgrastAnnotation}
\title{
Annotated the sequence on MGRAST
}
\description{
A completely Annotation process integrates login to MGRAST, data uploading and job submisstion.
}
\usage{
runMgrastAnnotation(login.info, fn, new_project)
}
\arguments{
  \item{login.info}{
web information while login to MGRAST , more details see \code{\link{loginMgrast}}  
}
  \item{fn}{
character,the sequence to analyze 
}
  \item{new_project}{
character, the project name
}
}
\details{
This function achieve a process contains login to MGRAST, data uploading and job submmition, more easier interface with MGRAST for sequence annotation. For more details of the each step, see \code{\link{loginMgrast}}, \code{\link{uploadMgrast}}, \code{\link{submitMgrastJob}}.
}
\value{
The MGRAST id for your job as well as the account.
}

\seealso{
\code{\link{loginMgrast}}, \code{\link{uploadMgrast}}, \code{\link{submitMgrastJob}}
}
\examples{
##
# runMgrastAnnotation(login.info, seqfiles,new_project)
}

