# Finding a Less Risky Alternative to the S&P 500

## Motivation
I wanted to find a portfolio of stocks that could deliver clients lower risk than the S&P 500 index while maintaining comparable returns.

## Tools
I used the programming language R to conduct this analysis, including:
- Calculating log returns, risk-adjusted returns (Sharpe ratios), and kurtosis.
- Minimizing portfolio variance to determine optimal weights for reducing volatility.

## Conducting Preliminary Study
- **Stock Selection**: Retrieved 25 stocks from Yahoo Finance with favorable technical indicators since the Great Recession of 2007-2009:
  - Stocks where the 50-week moving average stayed above the 200-week moving average for most of the chart's history.
  - Stocks that either bounced off or never touched the 200-week moving average, indicating long-term upward strength.
- **Annual Log Returns**: Calculated mean annual log returns for each stock (2009-2022).
- **Risk-Adjusted Returns**: Used the Sharpe ratio:

$$
S = \frac{R_p - R_f}{\sigma_p}
$$

  where:
  - $$\( R_p \)$$: Portfolio return
  - $$\( R_f \)$$: Risk-free rate (2.28%)
  - $$\( \sigma_p \)$$: Portfolio standard deviation

- **Excess Kurtosis**: Measured fatter tail distributions:

$$
K_{excess} = K - 3
$$

  where $$\( K \)$$ is kurtosis, though the sample size limited strong inferences.
- **Top 10 Stocks**: Selected based on the highest Sharpe ratios for better curation and asset management.

## Portfolio Return and Volatility
- Constructed an equal-weighted portfolio of 10 stocks:

$$
w_i = \frac{1}{10}
$$

  for each stock.

- **Performance**:
  - **Portfolio Return**:

$$
\mu_p = \sum_{i=1}^{n} w_i \mu_i
$$

  achieving 17.64% annualized log return, higher than the S&P 500 (10.51%).

  - **Portfolio Volatility**:

$$
\sigma_p = \sqrt{w^T \Sigma w}
$$

  achieving 0.1199, lower than any individual stock and slightly below the S&P 500 (0.1205).

## Optimizing the Portfolio
- **Objective**: Construct a portfolio with lower risk than the S&P 500 while maintaining similar returns.
- **Efficient Frontier**:
  - Used Modern Portfolio Theory to plot the convex envelope of possible portfolios by varying weights.
  - Constructed portfolios $$\( w \)$$ and $$\( v \)$$ to define the envelope, enabling interpolation of other portfolios via convex combinations:

$$
a w + (1-a) v
$$

  where $$\( a \)$$ is a constant.

- **Global Minimum Variance Portfolio (GMVP)**:
  - Minimized portfolio covariance:

$$
\sigma_p^2 = w^T \Sigma w
$$

  - Solved for weights $$\( w \)$$ using quadratic programming to achieve the least possible variance:

$$
w = \frac{\Sigma^{-1} \mathbf{1}}{\mathbf{1}^T \Sigma^{-1} \mathbf{1}}
$$

  - GMVP achieved a mean annual return of 9.79% with **drastically reduced volatility of 0.0415**.

## Conclusion
The analysis successfully identified a portfolio with:
- Lower risk than the S&P 500.
- Comparable returns through strategic stock selection and weight optimization.
- The results demonstrate the potential to significantly reduce portfolio risk while maintaining reasonable returns.
- This project highlights the effectiveness of data-driven optimization techniques in outperforming benchmarks like the S&P 500 on a risk-adjusted basis. 

Future improvements will focus on:
- Identifying stocks to achieve returns of 20-25% annually while maintaining variance around 0.04-0.05.
- Restricting short sales for easier portfolio management.

---

### Appendix

#### 25 Asset Candidates

- AAPL, MSFT, GOOGL, AMZN, BRK-B, NVDA, V, MA, WMT, HD, TMO, COST, ABT, ADBE, DHR, NKE, INTC, AMD, TMUS, AMGN, LMT, INTU, PGR, ICE, DG.


