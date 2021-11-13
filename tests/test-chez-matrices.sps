#!/usr/bin/env scheme-script
;; -*- mode: scheme; coding: utf-8 -*- !#
;; Copyright (c) 2021 Liam Packer
;; SPDX-License-Identifier: MIT
#!r6rs

(import (rnrs (6)) (srfi :64 testing) (chez-matrices))

(test-begin "simple matrix operations")
(let ((m (make-matrix 4 5 )))
  (test-equal "make-matrix" '#(#(0 0 0 0 0)
                               #(0 0 0 0 0)
                               #(0 0 0 0 0)
                               #(0 0 0 0 0)) m)
  (test-equal "num-dims" (matrix-num-dims m) 2)
  (test-equal "matrix-rows" (matrix-rows m) 4)
  (test-equal "matrix-cols" (matrix-cols m) 5)
  (matrix-set! m 3 2 10)
  (test-equal "matrix-ref" (matrix-ref m 3 2) 10)
  (test-equal "matrix-set" '#(#(0 0 0 0 0)
                              #(0 0 0 0 0)
                              #(0 0 0 0 0)
                              #(0 0 10 0 0)) m)
  (test-assert "matrix?" (and
                          (matrix? m)
                          (not (matrix? 0))
                          (not (matrix? '#(0 0)))
                          (matrix? '#(#(#(0))))))
  (test-assert "matrix-identity+matrix="
    (matrix= '#(#(1 0 0) #(0 1 0) #(0 0 1)) (matrix-identity 3)))
  (matrix-set-diagonal! m '#(3 3 3 3))
  (test-equal "matrix-set-diagonal"
    '#(#(3 0 0 0 0)
       #(0 3 0 0 0)
       #(0 0 3 0 0)
       #(0 0 10 3 0)) m)

  (matrix-set-row! m 0 '#(4 4 4 4 4))
  (test-equal "matrix-set-diagonal"
    '#(#(4 4 4 4 4)
       #(0 3 0 0 0)
       #(0 0 3 0 0)
       #(0 0 10 3 0)) m)

  (matrix-set-col! m 3 '#(5 5 5 5))
  (test-equal "matrix-set-diagonal"
    '#(#(4 4 4 5 4)
       #(0 3 0 5 0)
       #(0 0 3 5 0)
       #(0 0 10 5 0)) m)

  (test-equal "tr"
    (tr m) (+ 4 3 3 5))

  ;; (test-equal "det"
  ;;   (det) (* 4 3 3 5))

  )
(test-end)

;; (exit (if (zero? (test-runner-fail-count (test-runner-get))) 0 1))
