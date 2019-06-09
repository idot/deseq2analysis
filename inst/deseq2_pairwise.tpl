
## DESeq2 `r tit`
### p.adj < `r pcut` |lFC| > `r lfc` 
  
  * `r getPosExplanation(tit)`
  * `r getNegExplanation(tit)`
  
  
### MA plot
```{r ma-<% titfiles %>, fig.cap="MA Plot", include=TRUE}
   plotMAPlot(comp, pcut, lfc, comparisonFoldchange)
```


### Vulcano plot
```{r vulcano-<% titfiles %>, fig.cap="Vulcano Plot", include=TRUE}
   plotVulcano(comp, pcut, lfc, comparisonFoldchange)
```


### p-value distribution
```{r pvaldist-<% titfiles %>, fig.cap="p-Value distribution", include=TRUE}
   plotPvalDist(comp, pcut)
```


### independent filtering
```{r independent-<% titfiles %>, fig.cap="independent filtering", include=TRUE}
   plotIndependentFiltering(comp, pcut, lfc, filterThreshold)
```





