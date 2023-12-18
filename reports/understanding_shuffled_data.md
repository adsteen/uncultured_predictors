# Understanding shuffled data sets

It looks to me like the shuffled data set is weirdly similar to the
actual data set, in a way that makes me think I’ve misunderstood
something. First let’s load the unshuffled and shuffled data:

``` r
library(tidyverse)
theme_set(theme_classic())

d_real <- read_csv("../data/ashley_rerun_SMOTE/GEM-actual-bootstrapped-by-phylum-pathway-LASSO-stats-SMOTE-full.csv")
d_null <- read_csv("../data/ashley_rerun_SMOTE/GEM-shuffled-bootstrapped-by-phylum-pathway-LASSO-stats-SMOTE-full.csv")
```

There are clearly plenty of differences between the two data frames:

``` r
all.equal(d_real, d_null)
```

     [1] "Component \"phylum_name\": 1120 string mismatches"                              
     [2] "Component \"feature_name\": 1735 string mismatches"                             
     [3] "Component \"coef\": Mean relative difference: 1.305804"                         
     [4] "Component \"coef_sd\": Mean relative difference: 0.7924054"                     
     [5] "Component \"lower_95\": 'is.NA' value mismatch: 335 in current 590 in target"   
     [6] "Component \"upper_95\": 'is.NA' value mismatch: 335 in current 590 in target"   
     [7] "Component \"count\": Mean relative difference: 0.6968937"                       
     [8] "Component \"significant\": 474 element mismatches"                              
     [9] "Component \"permutation_importance\": Mean relative difference: 1.819572"       
    [10] "Component \"t-statistic\": 'is.NA' value mismatch: 335 in current 590 in target"
    [11] "Component \"p-value\": 'is.NA' value mismatch: 335 in current 590 in target"    

How are the data distributed?

``` r
coefs <- d_real |> 
  select(`p-value`) |> 
  mutate(source = "actual") |> 
  rbind(d_null |> select(`p-value`) |> mutate(source = "shuffled")) |> 
  filter(!is.na(`p-value`))
ggplot(coefs, aes(x=`p-value`, color=source)) +
  geom_density() 
```

![](understanding_shuffled_data_files/figure-commonmark/unnamed-chunk-3-1.png)

The p-values seem not identically distributed, which is good, but they
also seem very nearly identically distributed. Also, what’s up with the
hump in p-value density around p ~ 0.3?

Note that on a log scale the difference in p-value distribution is a lot
more evident:

``` r
ggplot(coefs, aes(x=log10(`p-value`), color=source)) +
  geom_density()
```

![](understanding_shuffled_data_files/figure-commonmark/unnamed-chunk-4-1.png)

# So what’s the problem?

The problem is that when I throw out insignificant rows, I end up
retaining exactly the same number of rows in the shuffled data set as in
the non-shuffled data set.

I’ve written some functions to do this filtering:

``` r
process_data <- function(d,
                         p.cutoff=0.01, 
                         adjust.p = TRUE, 
                         p.method = "holm"#,
                         #font.size = 10,
                         #legend.size = 10,
                         #extra.title.text = NULL
                         ) {
  
  d_filtered <- d |>  
    # adjust_p() replaces p-values of NA to 1, which is important for the holm correction
    adjust_p(d, p.cutoff, adjust.p, p.method = "holm") |> 
    filter(`p-value` <= p.cutoff) |> 
    # Note that LASSO likes to push coefficients to zero, which means p-value is NA
    mutate(towards.cultured = coef > 0)
  
  d_filtered
  }

adjust_p <- function(d, 
                     p.cutoff=0.01, 
                     p_adjust = FALSE, 
                     p.method = "holm",
                     font.size = 10,
                     legend.size = 10,
                     extra.title.text = NULL) {
  
  # Set rows where slope is 0 or p-value is NA to have a p-value of 1
  d <- d |> 
    mutate(`p-value` = case_when(coef == 0 ~ 1,
                                 is.na(`p-value`) ~ 1, # I think that is redundant with the previous line
                                 TRUE ~ `p-value`))
  
  # Adjust p values
  if(p_adjust) {
    d <- d |> 
      mutate(`p-value` = p.adjust(`p-value`, method = p.method))
  }
  d
}
```

## More “significant” associations in the shuffled data set than the real one

Let’s look at the data frames when I adjust p-values (using alpha=0.01)
and throw out insignificant rows:

``` r
d_real_filt <- process_data(d_real)
nrow(d_real_filt)
```

    [1] 635

``` r
d_null_filt <- process_data(d_null)
nrow(d_null_filt)
```

    [1] 658

I actually see more lines retained for the shuffled data than the null
data, although it is pretty close in both cases.

FWIW the ‘yield’ (i.e. number of rows retained) for the real data is
36.29%, while the yield for the shuffled data is 37.6.

If the statistical technique is working correctly and there are real
statistical associations, shouldn’t the shuffled data have far fewer
“significant” values than the unshuffled data frame?

FWIW this is not an artifact of p-value adjustment - you get the same
answer without it:

``` r
nrow(d_real |> filter(`p-value` < 0.01))
```

    [1] 813

``` r
nrow(d_null |> filter(`p-value` < 0.01))
```

    [1] 960
