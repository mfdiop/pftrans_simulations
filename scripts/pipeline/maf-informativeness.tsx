import React, { useState } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, ReferenceLine } from 'recharts';

const MAFInformativenessVisualization = () => {
  const [errorRate, setErrorRate] = useState(0.01);
  
  // Generate data for different MAF values
  const generateData = () => {
    const data = [];
    for (let maf = 0.001; maf <= 0.5; maf += 0.005) {
      const p = maf;
      
      // Informativeness: heterozygosity 2p(1-p)
      const informativeness = 2 * p * (1 - p);
      
      // P(match by chance) for homozygous
      const pMatchByChance = p * p + (1 - p) * (1 - p);
      
      // P(match if IBD)
      const pMatchIBD = 1 - errorRate;
      
      // LOD score for a match
      const lod = Math.log10(pMatchIBD / Math.max(pMatchByChance, 1e-10));
      
      // Reliability score (inverse of coefficient of variation)
      // For binomial: CV = sqrt((1-p)/np) where n is sample size
      // Assume n=100 samples
      const n = 100;
      const reliability = Math.sqrt(n * p * (1 - p)) / p;
      
      // Combined score: informativeness * reliability * LOD
      const combinedScore = informativeness * Math.min(reliability, 10) * lod / 10;
      
      data.push({
        maf: parseFloat(maf.toFixed(3)),
        informativeness: parseFloat(informativeness.toFixed(4)),
        lod: parseFloat(lod.toFixed(4)),
        reliability: parseFloat(Math.min(reliability, 10).toFixed(4)),
        combinedScore: parseFloat(Math.max(combinedScore, 0).toFixed(4))
      });
    }
    return data;
  };
  
  const data = generateData();
  
  return (
    <div className="w-full max-w-6xl mx-auto p-6 bg-gray-50">
      <div className="bg-white rounded-lg shadow-lg p-6 mb-6">
        <h2 className="text-2xl font-bold mb-2 text-gray-800">
          Why Filter on MAF? The Informativeness Trade-off
        </h2>
        <p className="text-gray-600 mb-4">
          Very rare alleles are theoretically informative but practically unreliable due to 
          genotyping errors and sampling variance.
        </p>
        
        <div className="mb-6">
          <label className="block text-sm font-medium text-gray-700 mb-2">
            Genotyping Error Rate: {(errorRate * 100).toFixed(1)}%
          </label>
          <input
            type="range"
            min="0.001"
            max="0.05"
            step="0.001"
            value={errorRate}
            onChange={(e) => setErrorRate(parseFloat(e.target.value))}
            className="w-full"
          />
        </div>
      </div>
      
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 mb-6">
        <div className="bg-white rounded-lg shadow-lg p-6">
          <h3 className="text-lg font-semibold mb-4 text-gray-800">
            Informativeness (Heterozygosity)
          </h3>
          <ResponsiveContainer width="100%" height={250}>
            <LineChart data={data}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis 
                dataKey="maf" 
                label={{ value: 'Minor Allele Frequency', position: 'insideBottom', offset: -5 }}
              />
              <YAxis label={{ value: 'Score', angle: -90, position: 'insideLeft' }} />
              <Tooltip />
              <Line type="monotone" dataKey="informativeness" stroke="#3b82f6" strokeWidth={2} dot={false} />
              <ReferenceLine x={0.01} stroke="red" strokeDasharray="3 3" label="Min MAF" />
            </LineChart>
          </ResponsiveContainer>
          <p className="text-sm text-gray-600 mt-2">
            Peak at MAF=0.5, but even moderate frequencies are informative
          </p>
        </div>
        
        <div className="bg-white rounded-lg shadow-lg p-6">
          <h3 className="text-lg font-semibold mb-4 text-gray-800">
            LOD Score (for match)
          </h3>
          <ResponsiveContainer width="100%" height={250}>
            <LineChart data={data}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis 
                dataKey="maf" 
                label={{ value: 'Minor Allele Frequency', position: 'insideBottom', offset: -5 }}
              />
              <YAxis label={{ value: 'LOD', angle: -90, position: 'insideLeft' }} />
              <Tooltip />
              <Line type="monotone" dataKey="lod" stroke="#10b981" strokeWidth={2} dot={false} />
              <ReferenceLine x={0.01} stroke="red" strokeDasharray="3 3" label="Min MAF" />
            </LineChart>
          </ResponsiveContainer>
          <p className="text-sm text-gray-600 mt-2">
            Rarer alleles give stronger LOD scores when they match
          </p>
        </div>
        
        <div className="bg-white rounded-lg shadow-lg p-6">
          <h3 className="text-lg font-semibold mb-4 text-gray-800">
            Reliability (sampling stability)
          </h3>
          <ResponsiveContainer width="100%" height={250}>
            <LineChart data={data}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis 
                dataKey="maf" 
                label={{ value: 'Minor Allele Frequency', position: 'insideBottom', offset: -5 }}
              />
              <YAxis label={{ value: 'Score', angle: -90, position: 'insideLeft' }} />
              <Tooltip />
              <Line type="monotone" dataKey="reliability" stroke="#f59e0b" strokeWidth={2} dot={false} />
              <ReferenceLine x={0.01} stroke="red" strokeDasharray="3 3" label="Min MAF" />
            </LineChart>
          </ResponsiveContainer>
          <p className="text-sm text-gray-600 mt-2">
            Very rare alleles have unstable frequency estimates
          </p>
        </div>
        
        <div className="bg-white rounded-lg shadow-lg p-6">
          <h3 className="text-lg font-semibold mb-4 text-gray-800">
            Combined Score (Overall Usefulness)
          </h3>
          <ResponsiveContainer width="100%" height={250}>
            <LineChart data={data}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis 
                dataKey="maf" 
                label={{ value: 'Minor Allele Frequency', position: 'insideBottom', offset: -5 }}
              />
              <YAxis label={{ value: 'Score', angle: -90, position: 'insideLeft' }} />
              <Tooltip />
              <Line type="monotone" dataKey="combinedScore" stroke="#8b5cf6" strokeWidth={2} dot={false} />
              <ReferenceLine x={0.01} stroke="red" strokeDasharray="3 3" label="Min MAF" />
              <ReferenceLine x={0.05} stroke="green" strokeDasharray="3 3" label="Optimal" />
            </LineChart>
          </ResponsiveContainer>
          <p className="text-sm text-gray-600 mt-2">
            <span className="font-semibold">Optimal range: MAF 0.05-0.30</span> - balances all factors
          </p>
        </div>
      </div>
      
      <div className="bg-white rounded-lg shadow-lg p-6">
        <h3 className="text-lg font-semibold mb-3 text-gray-800">Key Insights:</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <div className="border-l-4 border-red-500 pl-4">
            <h4 className="font-semibold text-red-700">Very Rare (MAF &lt; 0.01)</h4>
            <ul className="text-sm text-gray-600 mt-2 space-y-1">
              <li>• High LOD when matched</li>
              <li>• But: unstable frequency estimates</li>
              <li>• Dominated by genotyping errors</li>
              <li>• Risk of false positives/negatives</li>
            </ul>
          </div>
          
          <div className="border-l-4 border-green-500 pl-4">
            <h4 className="font-semibold text-green-700">Optimal (MAF 0.05-0.30)</h4>
            <ul className="text-sm text-gray-600 mt-2 space-y-1">
              <li>• Still informative for IBD</li>
              <li>• Reliable frequency estimates</li>
              <li>• Robust to genotyping errors</li>
              <li>• Best signal-to-noise ratio</li>
            </ul>
          </div>
          
          <div className="border-l-4 border-blue-500 pl-4">
            <h4 className="font-semibold text-blue-700">Common (MAF &gt; 0.4)</h4>
            <ul className="text-sm text-gray-600 mt-2 space-y-1">
              <li>• Very reliable estimates</li>
              <li>• But: matches often by chance</li>
              <li>• Low LOD scores</li>
              <li>• Still useful with many SNPs</li>
            </ul>
          </div>
          
          <div className="border-l-4 border-purple-500 pl-4">
            <h4 className="font-semibold text-purple-700">Practical Recommendation</h4>
            <ul className="text-sm text-gray-600 mt-2 space-y-1">
              <li>• Filter: 0.01 &lt; MAF &lt; 0.99</li>
              <li>• Weight by informativeness</li>
              <li>• Use LOD scores for segments</li>
              <li>• Quality over quantity!</li>
            </ul>
          </div>
        </div>
      </div>
      
      <div className="mt-6 bg-amber-50 border-l-4 border-amber-500 p-4 rounded">
        <h4 className="font-semibold text-amber-800 mb-2">⚠️ The Bottom Line:</h4>
        <p className="text-sm text-gray-700">
          We filter on MAF not because rare alleles aren't informative, but because 
          <span className="font-semibold"> very rare alleles are unreliable</span>. The optimal 
          strategy uses moderately rare to common alleles (MAF 0.01-0.5) where we can trust 
          the data, then weights them by informativeness to emphasize the rarer ones.
        </p>
      </div>
    </div>
  );
};

export default MAFInformativenessVisualization;
