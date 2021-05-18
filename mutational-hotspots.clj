(ns guertu.mutational-hotspots
  (:require [dk.ative.docjure.spreadsheet :as docjure]))


;; This is a modfied version of code for Cui et al. 2020: https://www.nature.com/articles/s41467-019-14099-w
;; Read in genomic data
;;
;; Results in a list of mutation locations
(def mutations-loc  (->> (docjure/load-workbook "SNP_matrix.xlsx")
                         (docjure/select-sheet "SNPs")
                         (docjure/select-columns {:A :id :CI :pos})
                         (map :pos)
                         (map int)
                         sort))


;; Length of the reference genome AE107334: https://www.ncbi.nlm.nih.gov/nuccore/AE017334
(def chromo-length 5227419) ;; at NCBI ranges from position 1 to 5227419.

;;
;; Function to generate a random list of mutations.
;;
(defn fake-genome [mutations-loc]
    (loop [mutations #{}]
                           (if (< (count mutations) (count mutations-loc))
                             (recur (conj mutations (rand-int chromo-length)))
                             (sort mutations))))

(defn potential-clusters [mutpos-coll chromo-length window]
  "This function can return clusters that are subsets of other clusters of a larger size N.
   Something to keep in mind that subsets of a larger cluster might also be a significant
   cluster (possibly more so than the larger cluster)."
  (let [w window
        clusters (reverse (sort-by count (set (map (fn [m] (filter #(<= m % (+ m w -1)) mutpos-coll)) mutpos-coll))))]
    clusters))

(defn find-kn [k n mutpos-coll chromo-length]
  "Lists all the clusters with K or more mutations within window size N"
  (sort-by first (map #(take k %) (filter #(<= k (count %)) (potential-clusters mutpos-coll chromo-length n)))))

mutations-loc
;;
;; Hotspot analysis
;;


;; STEPS:
;;
;; For K is 5 to 2,
;;
;; 1. Find all clusters of size K mutations in the observed data, in a window size of up to N = 2,000bp.
;;
;;
;; repeat until K = 2, then set window size to N = 200bp


;;
;; CLUSTER SIZE 5
;;
(println "Clusters of size k = 5 at a window size smaller or equal to 2,000 basepairs.")
(println (find-kn 5 2000 mutations-loc chromo-length)) ;
(println)

;;
;; CLUSTER SIZE 4
;;
(println "Clusters of size k = 4 at a window size smaller or equal to 2,000 basepairs.")
(println (find-kn 4 2000 mutations-loc chromo-length)) ;
(println)

;;
;; CLUSTER SIZE 3
;;
(println "Clusters of size k = 3 at a window size smaller or equal to 2,000 basepairs.")
(println (find-kn 3 2000 mutations-loc chromo-length)) ;
(println)

;;
;; CLUSTER SIZE 2 at 50
;;
(println "Clusters of size k = 2 at a window size smaller or equal to 50 basepairs.")
(println (find-kn 2 50 mutations-loc chromo-length)) ;
(println)



;;
;; SUMMARY:
;;
;; List the significantly clustered found clusters (with the gene and type of mutation listed)
;;
;;
;; [52546 52551 52613 52624 52641 52651 52652 52672] ;; rpoZ
;; [2326351 2326713 2328767 2329129 2329254]         ;; cutC, cutC, YPO2050, intergenic, YPO2051
;; [2577993 2577993]                                 ;; YPO2292-YPO2293
;; [2606316 2606316]                                 ;; 2606217..2606501	asr	YPO2318	-	acid shock protein
;; [358874 358875]                                   ;; aspA, aspA
;; [581519 581523]                                   ;; YPO0535-YPO0536
;; [3705090 3705104]                                 ;; YPO3321-YPO3322
