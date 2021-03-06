% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{yelp}
\alias{yelp}
\title{Yelp Dataset}
\format{\code{vocab} a vector of unique words in the corpus vocabulary.

\code{docs} a list of documents in the corpus. Each item (represents a
 document) is a matrix (2 X U) of word frequencies, where U represents the
 number of unique words in a document. Each column in the matrix represents
 a unique word in a document and contains
 \itemize{
   \item vocabulary-id. the index of the word in the vocabulary (starts with 0)
   \item frequency. the relative frequency of the word in the document
 }

\code{docs.metadata} a matrix of document (article) metadata, where each
row represents a document with
 \itemize{
   \item doc.id. a unique article id
   \item review.id.
   \item reviewer.id.
   \item rating. customer rating
   \item restaurant. restaurant name
   \item row.word.count. the number of words in the article
   \item category. the category of the review
 }

\code{cids} a vector of document collection ids

\code{class.labels} a vector of categories (classes) in the corpus

\code{collection.labels} a vector of collections in the corpus

\code{ds.name} the corpus name (string)

\code{num.docs} the number of documents in the corpus

\code{V} the vocabulary size}
\source{
Articles are downloaded from the
 \href{http://uilab.kaist.ac.kr/research/WSDM11}{link}
}
\usage{
\code{data("yelp")}
}
\description{
It is a collection of restaurant reviews from Yelp. This corpus consists of
24,310 reviews and 9,517 unique words.
}
\note{
Created on July 26, 2015
}
\seealso{
Other datasets: \code{\link{news}}, \code{\link{nips}}
}
\author{
Clint P. George
}
