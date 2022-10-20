---
    title: "FunctionPlots"
author: "Subhadip Pal"
date: "February 14, 2019"
output: html_document
runtime: shiny
---



```{r setup, include=FALSE} knitr::opts_chunk$set(echo = TRUE)
```



Let $f(y ; A)$ be the probability density function of the Modified Half Normal distribution then for any $m>0$,
$$ f(y ; A)\leq  \frac{m^{\frac{m}{ A+m}}   }{(A+m)}  {x^{\left(1 - \frac{m}{ A+m}\right)-1} \exp(-x)}$$

    $$ \frac{m}{(m+A)} \log(m) +\log\left(\Gamma\left(\frac{A}{(m+A)}\right)\right) -\log(A+m) $$
    $$ \frac{1}{A+m}+\frac{A}{(A+m)^2}\log(m) -\frac{A}{(A+m)^2}\psi\left(\frac{A}{(m+A)} \right) -\frac{1}{(m+A)} =   \frac{A}{(A+m)^2} \left[ \log(m) -\psi\left(\frac{A}{(m+A)} \right)\right] $$
    ## Inputs and Outputs

    You can embed Shiny inputs and outputs in your document. Outputs are automatically updated whenever inputs change.  This demonstrates how a standard R plot can be made interactive by wrapping it in the Shiny `renderPlot` function. The `selectInput` and `sliderInput` functions create the input widgets used to drive the plot.

```{r , echo=FALSE}
inputPanel(
    sliderInput("A_adjust", label = "The value of parameter A",
                min = 0.0002, max = 5, value = 1, step = 0.005),

    sliderInput("bw_adjust", label = "The value of parameter m",
                min = 0.00001, max = 5, value = 1, step = 0.0005),
    sliderInput("x_min", label = "The value of parameter x_min",
                min = 0.0002, max = .1, value = 1, step = 0.005),
    sliderInput("x_max", label = "The value of parameter x_max",
                min = .1, max = 5, value = 1, step = 0.005)
)

renderPlot({

    A=1;
    f_proposal<-function(x, A,m=NULL){

        #if(is.character(m)){m=(sqrt(gamma^2+8*(alpha-1)*beta)-gamma)/(4*beta)}
        if(is.null(m)){m=1}

        val=m^(m/(m+A)) * x^(-m/(m+A))* exp(-x)/(A+m)
        return(val)

    }


    f_unscaled<-function(x,A){
        if(sum(x)<0){return(0)}
        val=exp(-x)/((A+x))
        return(val)
    }


    eps=.00001

    m_opt<-function(A){
        g<-function(x){

            return(log(x)-digamma(A/(A+x)))
        }
        m_opt=uniroot(f = g, interval = c(eps, 10+2/A), tol = .00000000001 )
        return(m_opt$root)
    }


    m_op<-m_opt(input$A_adjust )







    library(ggplot2)
    p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
    p<-p + stat_function(fun = f_proposal,colour="orange",n=1000, args = list(A=input$A_adjust,  m=input$bw_adjust))
    p=p+stat_function(fun = f_unscaled,colour="blue",n=1000,args = list(A=input$A_adjust))  +xlim(input$x_min,input$x_max)
    p<-p+geom_vline(xintercept = m_op , col="purple") + geom_vline(xintercept = input$bw_adjust )
    p
})



```





```{r , echo=FALSE}
inputPanel(
    sliderInput("m_adjust", label = "The value of parameter m",
                min = 0.02, max = 2, value = 1, step = 0.05),
    sliderInput("A_adjust", label = "The value of parameter A",
                min = 0.002, max = 5, value = 1, step = 0.5)
)

renderPlot({

    A=2; m=1
    f_proposal<-function(x, A,m=NULL){

        #if(is.character(m)){m=(sqrt(gamma^2+8*(alpha-1)*beta)-gamma)/(4*beta)}
        if(is.null(m)){m=1}

        val=m^(m/(m+A)) * x^(-m/(m+A))* exp(-x)/(A+m)
        return(val)

    }


    f_unscaled<-function(x,A){
        if(sum(x)<0){return(0)}
        val=exp(-x)/((A+x))
        return(val)
    }


    library(ggplot2)

    p<-ggplot(NULL, aes(c(0,1))) +
        geom_area(stat = "function", fun = f_proposal, fill = "orange",color='black', xlim = c(0.00001, 1), alpha=.5,args = list( A=input$A_adjust, m=input$m_adjust) ) +
        geom_area(stat = "function", fun = f_unscaled, fill = "green",color='black', xlim = c(0.00001, 1), alpha=.4,args = list( A=input$A_adjust))
    #p<-p+geom_point(size=12, aes(x = input$m_adjust, y = f_unscaled(input$m_adjust,A = input$A_adjust)),col="#9933FF", alpha=.6)
    p
    #p<-p+geom_point(size=4, aes(x = input$m_adjust, y = f_unscaled(input$m_adjust,alpha = alpha,beta = beta,gamma = input$gamma_adjust)),col="#660066", alpha=.95)
    #p + annotate("segment", x = input$m_adjust, xend =  input$m_adjust, y = .95*f_unscaled( input$m_adjust,alpha = alpha,beta = beta,gamma = input$gamma_adjust) , yend = 0, colour = "#9933FF", size=2, alpha=0.6, arrow=arrow())




})

```

#The second type of inequalities

$$
    f(y\mid A)\propto  \frac{y^{-1}e^{-\frac{{A}}{y}}}{(1+y)} =   \frac{y^{-1}e^{-\frac{{A}}{y}}}{(1+m \left(\frac{y}{m} \right))}  \leq
\frac{y^{-1}e^{-\frac{{A}}{y}}}{(1+m) \left(\frac{y}{m} \right)^{\frac{m}{m+1}}}
=\frac{ m^{\frac{m}{m+1}}}{(1+m) } \left[y^{-\frac{m}{m+1}-1 }e^{-\frac{{A}}{y}} \right]
$$




    ```{r , echo=FALSE}
inputPanel(
    sliderInput("A_adjust1", label = "The value of parameter A",
                min = 0.0002, max = 20, value = 1, step = 0.005),

    sliderInput("bw_adjust1", label = "The value of parameter m",
                min = 0.00001, max = 20, value = 1, step = 0.005),
    sliderInput("x_min1", label = "The value of parameter x_min",
                min = 0.0002, max = .1, value = .1, step = 0.005),
    sliderInput("x_max1", label = "The value of parameter x_max",
                min = .1, max = 20, value = 1, step = 0.005)
)

renderPlot({

    A=1;
    f_proposal<-function(x, A,m=NULL){

        #if(is.character(m)){m=(sqrt(gamma^2+8*(alpha-1)*beta)-gamma)/(4*beta)}
        if(is.null(m)){m=1}

        val=m^(m/(m+1)) * x^(-1-m/(m+1))* exp(-A/x)/(1+m)
        return(val)

    }


    f_unscaled<-function(x,A){
        if(sum(x)<0){return(0)}
        val=exp(-A/x)/(x*(1+x))
        return(val)
    }

    library(ggplot2)
    p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
    p<-p + stat_function(fun = f_proposal,colour="orange",n=1000, args = list(A=input$A_adjust1,  m=input$bw_adjust1))
    p+stat_function(fun = f_unscaled,colour="blue",n=1000,args = list(A=input$A_adjust1))  +xlim(input$x_min1,input$x_max1)

})



```



#The Third type of inequalities
If we make a change of variable $y=\frac{\tau^2}{A}$ then,

$$
    f(y\mid A)\propto  \frac{y^{-1}e^{-\frac{{1}}{y}}}{(1+Ay)} =   \frac{y^{-1}e^{-\frac{{1}}{y}}}{(1+Am \left(\frac{y}{m} \right))}  \leq
\frac{y^{-1}e^{-\frac{{1}}{y}}}{(1+Am) \left(\frac{y}{m} \right)^{\frac{Am}{Am+1}}}
=\frac{ m^{\frac{Am}{Am+1}}}{(1+Am) } \left[y^{-\frac{Am}{Am+1}-1 }e^{-\frac{{1}}{y}} \right]
$$



    ```{r , echo=FALSE}
inputPanel(
    sliderInput("A_adjust2", label = "The value of parameter A",
                min = 0.0002, max = 20, value = 1, step = 0.0005),

    sliderInput("bw_adjust2", label = "The value of parameter m",
                min = 0.0002, max = 10, value = 1, step = 0.001),
    sliderInput("x_min2", label = "The value of parameter x_min",
                min = 0.0002, max = .1, value = .1, step = 0.005),
    sliderInput("x_max2", label = "The value of parameter x_max",
                min = .001, max = 20, value = 1, step = 0.005)
)

renderPlot({

    A=1;eps=.00001
    f_proposal<-function(x, A,m=NULL){

        #if(is.character(m)){m=(sqrt(gamma^2+8*(alpha-1)*beta)-gamma)/(4*beta)}
        if(is.null(m)){m=1}
        Am=A*m
        val=m^(Am/(Am+1)) * x^(-1-Am/(Am+1))* exp(-1/x)/(1+Am)
        return(val)
    }


    f_unscaled<-function(x,A){
        if(sum(x)<0){return(0)}
        val=exp(-1/x)/((x)*(1+A*x))
        return(val)
    }


    m_opt<-function(A){
        g<-function(x){
            Am=A*x
            return(log(x)+digamma(Am/(1+Am)))
        }
        m_opt=uniroot(f = g, interval = c(eps, 10+2/A) )
        return(m_opt$root)
    }


    f<-function(m){
        Am=A*m
        frac=Am/(1+Am)
        val=frac*log(m)+lgamma(frac)-log(Am+1)
        return(val)
    }

    m_op<-m_opt(input$A_adjust2)
    library(ggplot2)
    p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
    p<-p + stat_function(fun = f_proposal,colour="orange",n=1000, args = list(A=input$A_adjust2,  m=input$bw_adjust2))
    p<-p+stat_function(fun = f_unscaled,colour="blue",n=1000,args = list(A=input$A_adjust2))  +xlim(input$x_min2,input$x_max2)
    p<-p+geom_vline(xintercept = m_op , col="purple") + geom_vline(xintercept = input$bw_adjust2 )
    p

})









```

# fourth kind of ineq
Let $f(y ; A)$ be the probability density function of the Modified Half Normal distribution then for any $m>0$,
$$ f(y ; A)= \frac{e^{-Ay}}{1+y}\leq  \frac{m^{\frac{m}{ 1+m}}   }{(1+m)}  {y^{\left(1 - \frac{m}{ 1+m}\right)-1} \exp(-Ay)}$$

    $$ \frac{m}{(m+1)} \log(m) +\log\left(\Gamma\left(\frac{1}{(m+1)}\right)\right)- \frac{m}{(m+1)} \log(A) -\log(1+m) $$

    $$ \frac{1}{1+m}+\frac{1}{(1+m)^2}\log(m) -\frac{1}{(1+m)^2}\psi\left(\frac{1}{(m+1)} \right) -\frac{1}{(m+1)^2}\log(A)-\frac{1}{(m+1)} =   \frac{1}{(1+m)^2} \left[ \log(\frac{m}{A}) -\psi\left(\frac{1}{(m+1)} \right)\right] $$


    ## Inputs and Outputs

    You can embed Shiny inputs and outputs in your document. Outputs are automatically updated whenever inputs change.  This demonstrates how a standard R plot can be made interactive by wrapping it in the Shiny `renderPlot` function. The `selectInput` and `sliderInput` functions create the input widgets used to drive the plot.

```{r , echo=FALSE}
inputPanel(
    sliderInput("A_adjust4", label = "The value of parameter A",
                min = 0.0002, max = 5, value = 1, step = 0.005),

    sliderInput("bw_adjust4", label = "The value of parameter m",
                min = 0.00001, max = 5, value = 1, step = 0.0005),
    sliderInput("x_min4", label = "The value of parameter x_min",
                min = 0.0002, max = .1, value = 1, step = 0.005),
    sliderInput("x_max4", label = "The value of parameter x_max",
                min = .1, max = 5, value = 1, step = 0.005)
)

renderPlot({

    A=1;
    f_proposal<-function(x, A,m=NULL){

        #if(is.character(m)){m=(sqrt(gamma^2+8*(alpha-1)*beta)-gamma)/(4*beta)}
        if(is.null(m)){m=1}

        val=m^(m/(m+1)) * x^(-m/(m+1))* exp(-A*x)/(1+m)
        return(val)

    }


    m_opt<-function(A){
        g<-function(x){

            return(log(x/A)-digamma(1/(1+x)))
        }
        m_opt=uniroot(f = g, interval = c(eps, 10+2/A), tol = .00000000001 )
        return(m_opt$root)
    }

    eps=.000001
    f<-function(m){
        Am=A*m
        frac=1/(1+m)
        val=frac*log(m/A)+lgamma(frac)-log(m+1)
        return(val)
    }

    m_op<-m_opt(input$A_adjust4)





    f_unscaled<-function(x,A){
        if(sum(x)<0){return(0)}
        val=exp(-A*x)/((1+x))
        return(val)
    }

    library(ggplot2)
    p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
    p<-p + stat_function(fun = f_proposal,colour="orange",n=1000, args = list(A=input$A_adjust4,  m=input$bw_adjust4))
    p<-p+stat_function(fun = f_unscaled,colour="blue",n=1000,args = list(A=input$A_adjust4))  +xlim(input$x_min4,input$x_max4)
    p<-p+geom_vline(xintercept = m_op , col="purple") + geom_vline(xintercept = input$bw_adjust4 )
    p

})



```

