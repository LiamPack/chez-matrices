;; -*- mode: scheme; coding: utf-8 -*-
;; Copyright (c) 2021 Liam Packer
;; SPDX-License-Identifier: MIT
#!r6rs

(library (chez-matrices)
  (export make-matrix matrix-ref matrix-set! matrix-shape matrix-rank
          matrix-num-vals matrix-num-rows matrix-num-cols matrix? matrix-copy
          matrix-flatten vec->row-matrix matrix-ref-row matrix-set-row!
          matrix-set-diagonal! matrix-identity
          matrix-max matrix-min matrix-contract matrix-ref-col matrix-set-col!
          matrix-map mul T matrix-fold tr euclidean-norm cross-product
          theta phi linsolve matrix-inverse matrix= do-matrix matrix-fold
          matrix-generate)
  (import (rnrs (6))
          (only (srfi :133 vectors) vector-fold))

  (define (matrix-shape tensor)
    (let tensor-dims-rec ([tensor tensor] [dim-list (list)])
      (if (not (vector? (vector-ref tensor 0)))
          (reverse (cons (vector-length tensor) dim-list))
          (tensor-dims-rec (vector-ref tensor 0) (cons (vector-length tensor) dim-list)))))

  (define (matrix-rank tensor)
    (length (matrix-shape tensor)))

  (define-syntax do-matrix
    (syntax-rules ()
      [(_ m () rest ...) (begin rest ...)]
      [(_ m (() j ...) rest ...)
       (do ([i 0 (+ 1 i)])
           ((>= i 1))
         (do-matrix (matrix-ref m 0) (j ...) rest ...))]
      [(_ m (i j ...) rest ...)
       (do ([i 0 (+ 1 i)])
           ((>= i (car (matrix-shape m))))
         (do-matrix (matrix-ref m 0) (j ...) rest ...))]
      [(_ m ((i Upper) j ...) rest ...)
       (do ([i 0 (+ 1 i)])
           ((>= i Upper))
         (do-matrix (matrix-ref m 0) (j ...) rest ...))]
      [(_ m ((i Lower Upper) j ...) rest ...)
       (do ([i Lower (+ 1 i)])
           ((>= i Upper))
         (do-matrix (matrix-ref m 0) (j ...) rest ...))]
      [(_ m ret (i j ...) rest ...)
       (do ([i 0 (+ 1 i)])
           ((>= i (car (matrix-shape m))) ret)
         (do-matrix (matrix-ref m 0) (j ...) rest ...))]))

  (define-syntax make-matrix
    (syntax-rules ()
      [(_ i) (make-vector i)]
      [(_ i j)
       (do ([m (make-vector i)]
            [k 0 (+ k 1)])
           ((= k i) m)
         (vector-set! m k (make-vector j)))]
      [(_ i j k ...)
       (make-vector i (make-matrix j k ...))]))

  (define (matrix-ref m . ids)
    (if (= 1 (length ids))
        (vector-ref m (car ids))
        (apply matrix-ref (cons (vector-ref m (car ids)) (cdr ids)))))

  (define-syntax matrix-set!
    (syntax-rules ()
      [(_ m i val) (vector-set! m i val)]
      [(_ m i j k ... val) (matrix-set! (vector-ref m i) j k ... val)]))

  (define (matrix-num-vals m) (map * (matrix-shape m)))
  (define (matrix-num-rows a) (vector-length a))
  (define (matrix-num-cols a)
    (if (vector? (matrix-ref a 0))
        (vector-length (matrix-ref a 0))
        1))

  (define matrix?
    (lambda (x)
      (and (vector? x)
           (> (vector-length x) 0)
           (vector? (vector-ref x 0)))))

  ;; row matrix
  (define (vec->row-matrix v)
    (let ([w (make-vector 1)])
      (vector-set! w 0 v)
      w))

  (define (matrix-flatten m)
    (let* ([v (make-vector (* (matrix-num-cols m) (matrix-num-rows m)))])
      (do-matrix m v (i j)
                 (vector-set! v (+ j (* i (matrix-num-cols m)))
                              (matrix-ref m i j)))))

  (define (matrix-ref-row m i) (matrix-ref m i))
  (define (matrix-set-row! m i v) (matrix-set! m i v))

  ;; verry inefficient but clean!
  (define (matrix-ref-col m i)
    (matrix-ref-row (T m) i))
  (define (matrix-set-col! m i v)
    (do-matrix m (row ()) (matrix-set! m row i (vector-ref v row))))

  ;; assert len(v) = min(dimensions(m))
  (define (matrix-set-diagonal! m v)
    (let ([new-m (matrix-copy m)])
      (do ([i 0 (+ 1 i)])
          ((>= i (vector-length v)) m)
        (matrix-set! m i i (vector-ref v i)))))

  (define (matrix-identity dim)
    (let ([m (make-matrix dim dim)])
      (matrix-set-diagonal! m (make-vector dim 1))
      m))

  (define (matrix-contract m1 m2 i k)
    (do ([j 0 (+ 1 j)]
         [accum 0 (+ accum (* (matrix-ref m1 i j)
                              (matrix-ref m2 j k)))])
        ((>= j (matrix-num-cols m1)) accum)))

  (define matrix-map
    (lambda (f . As)
      (if (> (matrix-rank (car As)) 2)
          (apply vector-map (cons (lambda r (apply matrix-map (cons f r))) As))
          (apply vector-map (cons (lambda r (apply vector-map (cons f r))) As)))))
  
  (define (%matrix-fold1 f init m)
    (if (> (matrix-rank m) 1)
        (vector-fold f init (vector-map (lambda (x) (%matrix-fold1 f init x)) m))
        (vector-fold f init m)))

  (define (matrix-fold f init . ms)
    (fold-left f init (map (lambda (m) (%matrix-fold1 f init m)) ms)))

  (define-syntax matrix-generate
    (syntax-rules ()
      [(_  (lambda (i j ...) rest ...) n m ...)
       (let ([a (make-matrix n m ...)])
         (do-matrix a (i j ...)
                    (matrix-set! a i j ... ((lambda (i j ...) rest ...) i j ...)))
         a)]))

  ;; assert proper dims
  (define (matrix-copy a)
    (matrix-map (lambda (x) x) a))

  (define (matrix-copy! a b)
    (do-matrix a (i j) (matrix-set! a i j (matrix-ref b i j))))

  ;; throw error if dims of a1 and a2 are bad
  ;; case on a1/a2 being a constant or not
  (define (mul a1 a2)
    (cond
     [(and (matrix? a1) (matrix? a2))
      (let* ([m (make-matrix (matrix-num-rows a1) (matrix-num-cols a2))])
        (do-matrix m m (i k) (matrix-set! m i k (matrix-contract a1 a2 i k))))]
     [(and (matrix? a2) (number? a1))
      (matrix-map (lambda (x) (* a1 x)) a2)]
     [(and (matrix? a1) (number? a2))
      (matrix-map (lambda (x) (* a2 x)) a1)]))

  (define (tr m)
    (do ([i 0 (+ 1 i)]
         [accum 0 (+ accum (matrix-ref m i i))])
        ((>= i (min (matrix-num-cols m) (matrix-num-rows m))) accum)))

  (define (T a1)
    (let ([m (make-matrix (matrix-num-cols a1) (matrix-num-rows a1))])
      (do-matrix m m (i j)
                 (matrix-set! m j i (matrix-ref a1 i j)))))

  ;; assert len(v) = len(w) and both vectors and all that
  (define (dot v w)
    (matrix-ref (mul (T v) w) 0))

  (define matrix=
    (lambda (A B)
      (matrix-fold
       (lambda (x y) (and x y))
       #t
       (matrix-map = A B))))

  ;; Norms
  (define (matrix-max m)
    (matrix-fold max -inf.0 m))

  (define (matrix-min m)
    (matrix-fold min +inf.0 m))

  (define (frobenius-norm m)
    (matrix-fold
     +
     0 (matrix-map (lambda (x) (expt x 2)) m)))

  (define (euclidean-norm M)
    (cond
     [(matrix? M)
      (let ([ret 0])
        (matrix-map (lambda (x) (set! ret (+ ret (expt x 2)))) M))]
     [(vector? M)
      (do ([i 0 (+ 1 i)]
           [accum 0 (+ accum (expt (vector-ref M i) 2))])
          ((>= i (vector-length M)) (sqrt accum)))]))

  ;; Solve linear equation
  (define (linsolve input-A input-B)
    (when (and (matrix? input-A) (matrix? input-B))
      (let* ([A (matrix-copy input-A)]
             [B (matrix-copy input-B)]
             [nrows (matrix-num-rows A)]
             [ncols (matrix-num-cols A)])
        (let ([leading
               (lambda (col-idx)
                 (call/cc (lambda (break)
                            (do ([i col-idx (+ 1 i)])
                                ((or (>= col-idx nrows) (>= i ncols)) #f)
                              (when (not (zero? (matrix-ref A i col-idx)))
                                (break i))))))]
              [normalize-row
               (lambda (row col)
                 (let ([norm-val (matrix-ref A row col)])
                   (matrix-set-row! A row
                                    (vector-map (lambda (x) (/ x norm-val))
                                                (matrix-ref-row A row)))
                   (matrix-set! B row 0
                                (/ (matrix-ref B row 0) norm-val))))]
              [swap-rows
               (lambda (row-idx1 row-idx2)
                 (unless (= row-idx1 row-idx2)
                   (let ([tmp-row (matrix-ref-row A row-idx1)]
                         [tmp-b (matrix-ref B row-idx1 0)])
                     (matrix-set-row! A row-idx1 (matrix-ref-row A row-idx2))
                     (matrix-set-row! A row-idx2 tmp-row)
                     (matrix-set! B row-idx1 0 (matrix-ref B row-idx2 0))
                     (matrix-set! B row-idx2 0 tmp-b))))]
              [eliminate-row
               (lambda (row-idx col-idx basis-idx)
                 (let ([ratio (matrix-ref A row-idx col-idx)])
                   (do ([i col-idx (+ 1 i)])
                       ((>= i ncols))
                     (matrix-set! A row-idx i
                                  (- (matrix-ref A row-idx i)
                                     (* ratio (matrix-ref A basis-idx i)))))
                   (matrix-set! B row-idx 0
                                (- (matrix-ref B row-idx 0)
                                   (* ratio (matrix-ref B basis-idx 0))))))])
          (call/cc
           (lambda (break)
             (do ([j 0 (+ 1 j)])
                 ((>= j ncols) B)
               (let ([lead (leading j)])
                 (when (not lead)
                   (break #f))
                 (normalize-row lead j)
                 (do ([i 0 (+ 1 i)])
                     ((>= i nrows))
                   (unless (= i lead)
                     (eliminate-row i j lead)))
                 (swap-rows j lead)))))))))

  (define (matrix-inverse input-A)
    (when (and (= (matrix-num-cols input-A) (matrix-num-rows input-A)) (matrix? input-A))
      (let* ([A (matrix-copy input-A)]
             [nrows (matrix-num-rows A)]
             [ncols (matrix-num-cols A)]
             [inv-A (matrix-identity nrows)])
        (let ([leading
               (lambda (col-idx)
                 (call/cc
                  (lambda (break)
                    (do ([i col-idx (+ 1 i)])
                        ((or (>= col-idx nrows) (>= i ncols)) #f)
                      (when (not (zero? (matrix-ref A i col-idx)))
                        (break i))))))]
              [normalize-row
               (lambda (row col)
                 (let ([norm-val (matrix-ref A row col)])
                   (matrix-set-row! A row
                                    (vector-map (lambda (x) (/ x norm-val))
                                                (matrix-ref-row A row)))
                   (matrix-set-row! inv-A row
                                    (vector-map (lambda (x) (/ x norm-val))
                                                (matrix-ref-row inv-A row)))))]
              [swap-rows
               (lambda (row-idx1 row-idx2)
                 (unless (= row-idx1 row-idx2)
                   (let ([tmp-row (matrix-ref-row A row-idx1)]
                         [tmp-row-inv (matrix-ref-row inv-A row-idx1)])
                     (matrix-set-row! A row-idx1 (matrix-ref-row A row-idx2))
                     (matrix-set-row! A row-idx2 tmp-row)
                     (matrix-set-row! inv-A row-idx1 (matrix-ref-row inv-A row-idx2))
                     (matrix-set-row! inv-A row-idx2 tmp-row-inv))))]
              [eliminate-row
               (lambda (row-idx col-idx basis-idx)
                 (let ([ratio (matrix-ref A row-idx col-idx)])
                   (do ([i 0 (+ 1 i)])
                       ((>= i ncols))
                     (matrix-set! A row-idx i
                                  (- (matrix-ref A row-idx i)
                                     (* ratio (matrix-ref A basis-idx i))))
                     (matrix-set! inv-A row-idx i
                                  (- (matrix-ref inv-A row-idx i)
                                     (* ratio (matrix-ref inv-A basis-idx i)))))))])
          (call/cc
           (lambda (break)
             (do ([j 0 (+ 1 j)])
                 ((>= j ncols) inv-A)
               (let ([lead (leading j)])
                 (when (not lead)
                   (break #f))
                 (normalize-row lead j)
                 (do ([i 0 (+ 1 i)])
                     ((>= i nrows))
                   (unless (= i lead)
                     (eliminate-row i j lead)))
                 (swap-rows j lead)))
             ))))))

  ;; Using the physics notation of the azimuthal 2\pi angle being phi
  (define (phi v)
    (atan (vector-ref v 1)
          (vector-ref v 0)))

  (define (theta v)
    (acos (/ (vector-ref v 2)
             (euclidean-norm v))))

  ;; assert dim 3
  (define (cross-product v w)
    (let ([cross-vw (make-vector 3)])
      (apply
       (lambda (ax ay az bx by bz)
         (vector-set! cross-vw 0 (- (* ay bz)
                                    (* by az)))
         (vector-set! cross-vw 1 (- (- (* ax bz)
                                       (* bx az))))
         (vector-set! cross-vw 2 (- (* ax by)
                                    (* bx ay))))
       (append (vector->list v) (vector->list w)))
      cross-vw)))
