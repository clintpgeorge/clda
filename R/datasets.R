################################################################################
#
# Datasets: This file is part of clda
#
# Copyright (c) 2016  Clint P. George
#
# clda is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# clda is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

#' @name news
#'
#' @title 16 Newsgroups Dataset
#'
#' @description
#' It's a subset of the
#' \href{http://qwone.com/~jason/20Newsgroups}{20Newsgroups} dataset. This
#' corpus consists of 10,764 news articles and 9,208 unique words.
#'
#' @docType data
#'
#' @usage
#' \code{data("news")}
#'
#' @format
#' \code{vocab} a vector of unique words in the corpus vocabulary.
#'
#' \code{docs} a list of documents in the corpus. Each item (represents a
#'  document) is a matrix (2 X U) of word frequencies, where U represents the
#'  number of unique words in a document. Each column in the matrix represents
#'  a unique word in a document and contains
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)
#'    \item frequency. the relative frequency of the word in the document
#'  }
#'
#' \code{docs.metadata} a matrix of document (article) metadata, where each
#' row represents a document with
#'  \itemize{
#'    \item category. the category assigned to the article
#'    \item name. the name of the news article from the 20Newsgroups dataset
#'    \item doclength. the number of words in the article
#'    \item collection. the collection name of each article
#'  }
#'
#' \code{cids} a vector of document collection IDs
#'
#' \code{class.labels} a vector of categories (classes) in the corpus
#'
#' \code{collection.labels} a vector of collections in the corpus
#'
#' \code{ds.name} the corpus name (string)
#'
#' \code{num.docs} the number of documents in the corpus
#'
#' \code{V} the vocabulary size
#'
#'
#' @source Articles are downloaded via
#' \href{http://scikit-learn.org/}{scikit-learn}
#'
#'
#' @family datasets
#'
#' @note Created on July 26, 2015
#'
#' @author Clint P. George
NULL


#' @name nips
#'
#' @title NIPS 00-18 Dataset
#'
#' @description
#' It's a subset of the NIPS dataset (Chechik 2007). It consists of papers
#' published in proceedings 00 to 18 of the  Neural Information Processing
#' Systems (NIPS) conference (i.e. years from 1988 to 2005). This corpus
#' consists of 2,741 news articles and 9,156 unique words.
#'
#' @docType data
#'
#' @usage
#' \code{data("nips")}
#'
#' @format
#' \code{vocab} a vector of unique words in the corpus vocabulary.
#'
#' \code{docs} a list of documents in the corpus. Each item (represents a
#'  document) is a matrix (2 X U) of word frequencies, where U represents the
#'  number of unique words in a document. Each column in the matrix represents
#'  a unique word in a document and contains
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)
#'    \item frequency. the relative frequency of the word in the document
#'  }
#'
#' \code{docs.metadata} a matrix of document (article) metadata, where each
#' row represents a document with
#'  \itemize{
#'    \item doc.id. a unique article id
#'    \item file.name. the name of the article
#'    \item row.word.count. the number of words in the article
#'    \item collection. the collection name of each article
#'  }
#'
#' \code{cids} a vector of document collection IDs
#'
#' \code{class.labels} a vector of categories (classes) in the corpus
#'
#' \code{collection.labels} a vector of collections in the corpus
#'
#' \code{ds.name} the corpus name (string)
#'
#' \code{num.docs} the number of documents in the corpus
#'
#' \code{V} the vocabulary size
#'
#'
#' @source Articles are downloaded from the
#'  \href{http://robotics.stanford.edu/~gal/data.html}{link}
#'
#'
#' @family datasets
#'
#' @note Created on July 26, 2015
#'
#' @author Clint P. George
NULL


#' @name yelp
#'
#' @title Yelp Dataset
#'
#' @description
#' It is a collection of restaurant reviews from Yelp. This corpus consists of
#' 24,310 reviews and 9,517 unique words.
#'
#' @docType data
#'
#' @usage
#' \code{data("yelp")}
#'
#' @format
#' \code{vocab} a vector of unique words in the corpus vocabulary.
#'
#' \code{docs} a list of documents in the corpus. Each item (represents a
#'  document) is a matrix (2 X U) of word frequencies, where U represents the
#'  number of unique words in a document. Each column in the matrix represents
#'  a unique word in a document and contains
#'  \itemize{
#'    \item vocabulary-id. the index of the word in the vocabulary (starts with 0)
#'    \item frequency. the relative frequency of the word in the document
#'  }
#'
#' \code{docs.metadata} a matrix of document (article) metadata, where each
#' row represents a document with
#'  \itemize{
#'    \item doc.id. a unique article id
#'    \item review.id.
#'    \item reviewer.id.
#'    \item rating. customer rating
#'    \item restaurant. restaurant name
#'    \item row.word.count. the number of words in the article
#'    \item category. the category of the review
#'  }
#'
#' \code{cids} a vector of document collection ids
#'
#' \code{class.labels} a vector of categories (classes) in the corpus
#'
#' \code{collection.labels} a vector of collections in the corpus
#'
#' \code{ds.name} the corpus name (string)
#'
#' \code{num.docs} the number of documents in the corpus
#'
#' \code{V} the vocabulary size
#'
#'
#' @source Articles are downloaded from the
#'  \href{http://uilab.kaist.ac.kr/research/WSDM11}{link}
#'
#'
#' @family datasets
#'
#' @note Created on July 26, 2015
#'
#' @author Clint P. George
NULL

