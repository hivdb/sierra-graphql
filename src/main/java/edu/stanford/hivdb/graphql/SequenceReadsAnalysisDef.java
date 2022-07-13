/*

    Copyright (C) 2017 Stanford HIVDB team

    Sierra is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sierra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.stanford.hivdb.graphql;

import graphql.GraphQLException;
import graphql.schema.*;
import static graphql.Scalars.*;
import static graphql.schema.GraphQLObjectType.newObject;
import static graphql.schema.GraphQLInputObjectType.newInputObject;
import static graphql.schema.GraphQLCodeRegistry.newCodeRegistry;
import static graphql.schema.FieldCoordinates.coordinates;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;

import edu.stanford.hivdb.drugresistance.GeneDR;
import edu.stanford.hivdb.genotypes.BoundGenotype;
import edu.stanford.hivdb.genotypes.GenotypeResult;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.PositionCodonReads;
import edu.stanford.hivdb.seqreads.CutoffCalculator.CutoffKeyPoint;
import edu.stanford.hivdb.seqreads.GeneSequenceReads;
import edu.stanford.hivdb.seqreads.OneCodonReadsCoverage;
import edu.stanford.hivdb.seqreads.SequenceReads;
import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.utilities.SimpleMemoizer;
import edu.stanford.hivdb.utilities.ValidationResult;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.UntranslatedRegion;
import edu.stanford.hivdb.viruses.Virus;

import static edu.stanford.hivdb.graphql.MutationSetDef.*;
import static edu.stanford.hivdb.graphql.GeneDef.*;
import static edu.stanford.hivdb.graphql.StrainDef.*;
import static edu.stanford.hivdb.graphql.ValidationResultDef.*;
import static edu.stanford.hivdb.graphql.GeneSequenceReadsDef.*;
import static edu.stanford.hivdb.graphql.DrugResistanceDef.*;
import static edu.stanford.hivdb.graphql.SubtypeV2Def.*;
import static edu.stanford.hivdb.graphql.PositionCodonReadsDef.*;
import static edu.stanford.hivdb.graphql.MutationPrevalenceDef.*;
import static edu.stanford.hivdb.graphql.AlgorithmComparisonDef.*;
import static edu.stanford.hivdb.graphql.SequenceReadsHistogramDef.*;
import static edu.stanford.hivdb.graphql.SequenceReadsHistogramByCodonReadsDef.*;
import static edu.stanford.hivdb.graphql.DrugResistanceAlgorithmDef.*;
import static edu.stanford.hivdb.graphql.DescriptiveStatisticsDef.*;

public class SequenceReadsAnalysisDef {
	
	private static <VirusT extends Virus<VirusT>> DataFetcher<List<BoundGenotype<VirusT>>> makeSubtypesDataFetcher(VirusT virusIns) {
		return env -> {
			int first = env.getArgument("first");
			SequenceReads<VirusT> seqReads = env.getSource();
			GenotypeResult<VirusT> subtypeResult = seqReads.getSubtypeResult();
			if (subtypeResult == null) {
				return Collections.emptyList();
			}
			List<BoundGenotype<VirusT>> matches = subtypeResult.getAllMatches();
			first = Math.min(first, matches.size());
			return matches.subList(0, first);
		};
	};

	private static <VirusT extends Virus<VirusT>> DataFetcher<List<ValidationResult>> makeValidationResultsDataFetcher(VirusT virusIns) {
		return env -> {
			SequenceReads<VirusT> seqReads = env.getSource();
			Collection<String> includeGenes = env.getArgument("includeGenes");
			return seqReads.getValidationResults(Sets.newLinkedHashSet(includeGenes));
		};
	}

	private static <VirusT extends Virus<VirusT>> DataFetcher<List<GeneSequenceReads<VirusT>>> makeAllGeneSequenceReadsDataFetcher(VirusT virusIns) {
		return env -> {
			SequenceReads<VirusT> seqReads = env.getSource();
			Collection<String> includeGenes = env.getArgument("includeGenes");
			return seqReads.getAllGeneSequenceReads(Sets.newLinkedHashSet(includeGenes));
		};
	}

	private static <VirusT extends Virus<VirusT>> DataFetcher<List<GeneDR<VirusT>>> makeDrugResistanceDataFetcher(VirusT virusIns) {
		return env -> {
			SequenceReads<VirusT> seqReads = env.getSource();
			String algName = env.getArgument("algorithm");
			Collection<String> includeGenes = env.getArgument("includeGenes");
			List<GeneSequenceReads<VirusT>> allGeneSeqReads = seqReads.getAllGeneSequenceReads(Sets.newLinkedHashSet(includeGenes));
			return new ArrayList<>(GeneDR.newFromGeneSequenceReads(
				allGeneSeqReads,
				virusIns.getDrugResistAlgorithm(algName)
			).values());
		};
	};

	private static <VirusT extends Virus<VirusT>> DataFetcher<List<OneCodonReadsCoverage<VirusT>>> makeCodonReadsCoverageDataFetcher(VirusT virusIns) {
		return env -> {
			SequenceReads<VirusT> sr = env.getSource();
			Collection<String> includeGenes = env.getArgument("includeGenes");
			return sr.getCodonReadsCoverage(includeGenes);
		};
	};
	

	private static DataFetcher<String> internalJsonCodonReadsCoverageDataFetcher = env -> {
		SequenceReads<?> sr = env.getSource();
		Collection<String> includeGenes = env.getArgument("includeGenes");
		return Json.dumpsUgly(
			sr
			.getCodonReadsCoverage(includeGenes)
			.stream()
			.map(rc -> rc.extMap())
			.collect(Collectors.toList()));
	};
	
	
	public static <VirusT extends Virus<VirusT>> SequenceReads<VirusT> toSequenceReadsList(Map<String, Object> input) {
		String name = (String) input.get("name");
		if (name == null) {
			throw new GraphQLException("`name` is a required field but doesn't have value");
		}
		@SuppressWarnings("unchecked")
		Strain<VirusT> strain = (Strain<VirusT>) input.get("strain");
		if (strain == null) {
		 throw new GraphQLException("`strain` is a required field but doesn't have value");
		}
		List<PositionCodonReads<VirusT>> allReads = (
			((List<?>) input.get("allReads"))
			.stream()
			.map(pcr -> toPositionCodonReads(strain, (Map<?, ?>) pcr))
			.collect(Collectors.toList()));
		if (allReads == null) {
			throw new GraphQLException("`allReads` is a required field but doesn't have value");
		}
		
		List<?> inputUntransRegions = (List<?>) input.get("untranslatedRegions");
		List<UntranslatedRegion> untransRegions = (
			inputUntransRegions == null ? Collections.emptyList() : inputUntransRegions
			.stream()
			.map(utr -> new UntranslatedRegion(
				(String) (((Map<?, ?>) utr).get("name")),
				(Long) (((Map<?, ?>) utr).get("refStart")),
				(Long) (((Map<?, ?>) utr).get("refEnd")),
				(String) (((Map<?, ?>) utr).get("consensus"))
			))
			.collect(Collectors.toList())
		);
		
		return SequenceReads.fromCodonReadsTable(
			(String) input.get("name"),
			strain,
			allReads,
			untransRegions,
			(Double) input.get("maxMixtureRate"),
			(Double) input.get("minPrevalence"),
			(Long) input.get("minCodonReads"),
			(Long) input.get("minPositionReads")
		);
	}

	public static GraphQLInputObjectType iUntranslatedRegion = (
		newInputObject()
		.name("UntranslatedRegionInput")
		.description("Optional untranslated region data.")
		.field(field -> field
			.type(GraphQLString)
			.name("name")
			.description("Name of untranslated region."))
		.field(field -> field
			.type(GraphQLLong)
			.name("refStart")
			.description("Absolute position (1-based) where this untranslated region started."))
		.field(field -> field
			.type(GraphQLLong)
			.name("refEnd")
			.description("Absolute position (1-based) where this untranslated region ended."))
		.field(field -> field
			.type(GraphQLString)
			.name("consensus")
			.description("NA Consensus of this untranslated region."))
		.build()
	);

	public static SimpleMemoizer<GraphQLInputType> iSequenceReads = new SimpleMemoizer<>(
		name -> (
			newInputObject()
			.name("SequenceReadsInput")
			.field(field -> field
				.type(GraphQLString)
				.name("name")
				.description("An identifiable name for identifying the result from the returning list."))
			.field(field -> field
				.type(enumStrain.get(name))
				.name("strain")
				.description("Strain of this sequence."))
			.field(field -> field
				.type(new GraphQLList(iPositionCodonReads.get(name)))
				.name("allReads")
				.description("List of all reads belong to this sequence."))
			.field(field -> field
				.type(new GraphQLList(iUntranslatedRegion))
				.defaultValue(null)
				.name("untranslatedRegions")
				.description("Optional consensus information of untranslated regions to be added in assembled consensus sequence."))
			.field(field -> field
				.type(GraphQLFloat)
				.name("maxMixtureRate")
				.defaultValue(1.)
				.description(
					"The maximum allowed mixture percentage cutoff. " +
					"Default to one if this field was left empty or had a " +
					"negative number specified. Valid value ranges from 0 to 1."))
			.field(field -> field
				.type(GraphQLFloat)
				.name("minPrevalence")
				.defaultValue(0.)
				.description(
					"The minimal prevalence cutoff to apply on each **codon**. " +
					"Default to zero if this field was left empty or had a " +
					"negative number specified. Valid value ranges from 0 to 1."))
			.field(field -> field
				.type(GraphQLLong)
				.name("minCodonReads")
				.defaultValue(1L)
				.description(
					"The minimal read depth for **codons**. " +
					"Default to one if this field was left empty or had a " +
					"negative number specified."))
			.field(field -> field
				.type(GraphQLLong)
				.name("minPositionReads")
				.defaultValue(1L)
				.description(
					"The minimal read depth for **positions**. " +
					"Default to one if this field was left empty or had a " +
					"negative number specified."))
			.build()
		)
	);

	private static DataFetcher<Boolean> isTrimmedDataFetcher = env -> {
		OneCodonReadsCoverage<?> ocrc = env.getSource();
		return ocrc.isTrimmed();
	};
	
	private static GraphQLCodeRegistry oneCodonReadsCoverageCodeRegistry = newCodeRegistry()
		.dataFetcher(
			coordinates("OneCodonReadsCoverage", "isTrimmed"),
			isTrimmedDataFetcher
		)
		.build();

	public static SimpleMemoizer<GraphQLObjectType> oOneCodonReadsCoverage = new SimpleMemoizer<>(
		name -> (
			newObject()
			.name("OneCodonReadsCoverage")
			.field(field -> field
				.type(oGene.get(name))
				.name("gene")
				.description("Gene of this record.")
			)
			.field(field -> field
				.type(GraphQLLong)
				.name("position")
				.description("Codon position in this gene.")
			)
			.field(field -> field
				.type(GraphQLLong)
				.name("totalReads")
				.description("Total reads of this position.")
			)
			.field(field -> field
				.type(GraphQLBoolean)
				.name("isTrimmed")
				.description("This position is trimmed or not.")
			)
			.build()
		)
	);
	
	private static DataFetcher<Boolean> isAboveMixtureRateThresholdFetcher = env -> {
		CutoffKeyPoint cutoff = env.getSource();
		return cutoff.isAboveMixtureRateThreshold();
	};
	
	private static DataFetcher<Boolean> isBelowMinPrevalenceThresholdFetcher = env -> {
		CutoffKeyPoint cutoff = env.getSource();
		return cutoff.isBelowMinPrevalenceThreshold();
	};
	
	private static GraphQLCodeRegistry cutoffKeyPointCodRegistry = newCodeRegistry()
		.dataFetcher(
			coordinates("CutoffKeyPoint", "isAboveMixtureRateThreshold"),
			isAboveMixtureRateThresholdFetcher
		)
		.dataFetcher(
			coordinates("CutoffKeyPoint", "isBelowMinPrevalenceThreshold"),
			isBelowMinPrevalenceThresholdFetcher
		)
		.build();

	private static GraphQLObjectType oCutoffKeyPoint = newObject()
		.name("CutoffKeyPoint")
		.description("Related `mixtureRate` and `minPrevalence`.")
		.field(field -> field
			.name("mixtureRate")
			.type(GraphQLFloat)
			.description("Mixture rate of all included NAs")
		)
		.field(field -> field
			.name("minPrevalence")
			.type(GraphQLFloat)
			.description("Minimal prevalence of all included codons")
		)
		.field(field -> field
			.name("isAboveMixtureRateThreshold")
			.type(GraphQLBoolean)
			.description("Is current `mixtureRate` above the input threshold `maxMixtureRate`?")
		)
		.field(field -> field
			.name("isBelowMinPrevalenceThreshold")
			.type(GraphQLBoolean)
			.description("Is current `minPrevalence` below the input argument `minPrevalence`?")
		)
		.build();

	private static DataFetcher<List<Map<String, Object>>> boundMutPrevListDataFetcher = env -> {
		SequenceReads<?> seqReads = env.getSource();
		MutationSet<?> mutations = seqReads.getMutations();
		Collection<String> includeGenes = env.getArgument("includeGenes");
		return getBoundMutationPrevalenceList(mutations, Sets.newHashSet(includeGenes));
	};

	private static <VirusT extends Virus<VirusT>> DataFetcher<List<Map<String, Object>>> makeMutAlgCmpDataFetcher(VirusT virusIns) {
		return env -> {
			List<String> asiAlgs = env.getArgument("algorithms");
			List<Map<String, String>> customAlgs = env.getArgument("customAlgorithms");
			if (asiAlgs == null) { asiAlgs = Collections.emptyList(); }
			if (customAlgs == null) { customAlgs = Collections.emptyList(); }
			if (asiAlgs.isEmpty() && customAlgs.isEmpty()) {
				return Collections.emptyList();
			}
			asiAlgs = asiAlgs
				.stream().filter(alg -> alg != null)
				.collect(Collectors.toList());
			Map<String, String> customAlgs2 = customAlgs
				.stream()
				.filter(map -> map != null)
				.collect(Collectors.toMap(
					map -> map.get("name"),
					map -> map.get("xml"),
					(x1, x2) -> x2,
					LinkedHashMap::new
				));
			SequenceReads<VirusT> seqReads = env.getSource();
			return fetchAlgorithmComparisonData(virusIns, seqReads.getMutations(), asiAlgs, customAlgs2);
		};
	};
	
	public static <VirusT extends Virus<VirusT>> GraphQLCodeRegistry makeSequenceReadsCodeRegistry(VirusT virusIns) {
		return (
			newCodeRegistry()
			.dataFetcher(
				coordinates("SequenceReadsAnalysis", "subtypes"),
				makeSubtypesDataFetcher(virusIns)
			)
			.dataFetcher(
				coordinates("SequenceReadsAnalysis", "validationResults"),
				makeValidationResultsDataFetcher(virusIns)
			)
			.dataFetcher(
				coordinates("SequenceReadsAnalysis", "allGeneSequenceReads"),
				makeAllGeneSequenceReadsDataFetcher(virusIns)
			)
			.dataFetcher(
				coordinates("SequenceReadsAnalysis", "drugResistance"),
				makeDrugResistanceDataFetcher(virusIns)
			)
			.dataFetcher(
				coordinates("SequenceReadsAnalysis", "codonReadsCoverage"),
				makeCodonReadsCoverageDataFetcher(virusIns)
			)
			.dataFetcher(
				coordinates("SequenceReadsAnalysis", "internalJsonCodonReadsCoverage"),
				internalJsonCodonReadsCoverageDataFetcher
			)
			.dataFetcher(
				coordinates("SequenceReadsAnalysis", "histogram"),
				seqReadsHistogramDataFetcher
			)
			.dataFetcher(
				coordinates("SequenceReadsAnalysis", "histogramByCodonReads"),
				seqReadsHistogramByCodonReadsDataFetcher
			)
			.dataFetcher(
				coordinates("SequenceReadsAnalysis", "mutations"),
				newMutationSetDataFetcher(virusIns, "mutations")
			)
			.dataFetcher(
				coordinates("SequenceReadsAnalysis", "mutationPrevalences"),
				boundMutPrevListDataFetcher
			)
			.dataFetcher(
				coordinates("SequenceReadsAnalysis", "algorithmComparison"),
				makeMutAlgCmpDataFetcher(virusIns)
			)
			.dataFetchers(oneCodonReadsCoverageCodeRegistry)
			.dataFetchers(cutoffKeyPointCodRegistry)
			.dataFetchers(makeGeneSequenceReadsCodeRegistry(virusIns))
			.build()
		);
	}
			
	public static SimpleMemoizer<GraphQLObjectType> oSequenceReadsAnalysis = new SimpleMemoizer<>(
		virusName -> {
			GraphQLObjectType.Builder builder = newObject()
				.name("SequenceReadsAnalysis")
				.field(field -> field
					.type(GraphQLString)
					.name("name")
					.description("Name of this sequence."))
				.field(field -> field
					.type(oStrain)
					.name("strain")
					.description("Strain of this sequence."))
				.field(field -> field
					.type(GraphQLFloat)
					.name("cutoffSuggestionLooserLimit")
					.description(
						"Algorithm suggested minimal prevalence cutoff. " +
						"This cutoff is looser and may include more problematic mutations."))
				.field(field -> field
					.type(GraphQLFloat)
					.name("cutoffSuggestionStricterLimit")
					.description(
						"Algorithm suggested minimal prevalence cutoff. " +
						"This cutoff is stricter and include less problematic mutations."))
				.field(field -> field
					.type(new GraphQLList(oValidationResult))
					.name("validationResults")
					.argument(arg -> arg
						.type(new GraphQLList(GeneDef.enumGene.get(virusName)))
						.name("includeGenes")
						.defaultValue(Virus.getInstance(virusName).getDefaultIncludedGenes())
						.description("Genes to be included in the results")
					)
					.description("Validation results for the sequence reads."))
				.field(field -> field
					.type(GraphQLFloat)
					.name("actualMinPrevalence")
					.description(
						"The actual minimal prevalence cutoff applied on this sequence."
					))
				.field(field -> field
					.type(GraphQLFloat)
					.name("minPrevalence")
					.description(
						"The input minimal prevalence cutoff applied on this sequence."
					))
				.field(field -> field
					.type(GraphQLLong)
					.name("minCodonReads")
					.description(
						"The minimal codon count cutoff applied on this sequence."
					))
				.field(field -> field
					.type(GraphQLLong)
					.name("minPositionReads")
					.description(
						"The minimal read depth for each position of the sequence reads."
					))
				.field(field -> field
					.type(new GraphQLList(oGene.get(virusName)))
					.name("availableGenes")
					.description("Available genes found in the sequence reads."))
				.field(field -> field
					.type(new GraphQLList(oGeneSequenceReads.get(virusName)))
					.name("allGeneSequenceReads")
					.argument(arg -> arg
						.type(new GraphQLList(GeneDef.enumGene.get(virusName)))
						.name("includeGenes")
						.defaultValue(Virus.getInstance(virusName).getDefaultIncludedGenes())
						.description("Genes to be included in the results")
					)
					.description("List of sequence reads distinguished by genes."))
				.field(field -> field
					.type(new GraphQLList(oBoundSubtypeV2))
					.name("subtypes")
					.argument(arg -> arg
						.type(GraphQLInt)
						.name("first")
						.defaultValue(2)
						.description(
							"Fetch only the first nth closest subtypes. Default to 2."))
					.description(
						"List of HIV1 groups or subtypes, or HIV species. " +
						"Sorted by the similarity from most to least."))
				.field(field -> field
					.type(oBoundSubtypeV2)
					.name("bestMatchingSubtype")
					.description(
						"The best matching subtype."))
				.field(field -> field
					.type(GraphQLFloat)
					.name("maxMixtureRate")
					.description(
						"Maximum allowed mixture percentage specified by input."
					))
				.field(field -> field
					.type(GraphQLFloat)
					.name("mixtureRate")
					.description(
						"Post-filter nucleotide mixture percentage."))
				.field(field -> newMutationSet(virusName, field, "mutations", /* enableIncludedGenes= */true)
					.description("All mutations found in the sequence reads."))
				.field(field -> field
					.type(GraphQLInt)
					.name("mutationCount")
					.description("Number of mutations without counting unsequenced regions and multiple continuous deletions"))
				.field(field -> field
					.type(GraphQLInt)
					.name("unusualMutationCount")
					.description("Number of unusual mutations without counting unsequenced regions and multiple continuous deletions"))
				.field(field -> field
					.type(new GraphQLList(oDrugResistance.get(virusName)))
					.name("drugResistance")
					.argument(arg -> arg
						.name("algorithm")
						.type(oASIAlgorithm.get(virusName))
						.defaultValue(Virus.getInstance(virusName).getDefaultDrugResistAlgorithm().getName())
						.description("One of the built-in ASI algorithms."))
					.argument(arg -> arg
						.type(new GraphQLList(GeneDef.enumGene.get(virusName)))
						.name("includeGenes")
						.defaultValue(Virus.getInstance(virusName).getDefaultIncludedGenes())
						.description("Genes to be included in the results")
					)
					.description("List of drug resistance results by genes."))
				.field(oSeqReadsHistogramBuilder)
				.field(oSeqReadsHistogramByCodonReadsBuilder)
				.field(field -> field
					.name("readDepthStats")
					.type(oDescriptiveStatistics)
					.description("Descriptive statistics of read depth for all positions.")
				)
				.field(field -> field
					.name("readDepthStatsDRP")
					.type(oDescriptiveStatistics)
					.description("Descriptive statistics of read depth for drug resistance positions.")
				)
				.field(field -> field
					.name("codonReadsCoverage")
					.type(new GraphQLList(oOneCodonReadsCoverage.get(virusName)))
					.argument(arg -> arg
						.type(new GraphQLList(GeneDef.enumGene.get(virusName)))
						.name("includeGenes")
						.defaultValue(Virus.getInstance(virusName).getDefaultIncludedGenes())
						.description("Genes to be included in the results")
					)
					.description("Codon reads coverage.")
				)
				.field(field -> field
					.type(GraphQLString)
					.name("internalJsonCodonReadsCoverage")
					.argument(arg -> arg
						.type(new GraphQLList(GeneDef.enumGene.get(virusName)))
						.name("includeGenes")
						.defaultValue(Virus.getInstance(virusName).getDefaultIncludedGenes())
						.description("Genes to be included in the results")
					)
					.description(
						"Position codon reads in this gene sequence (json formated)."))
				.field(field -> field
					.type(new GraphQLList(oCutoffKeyPoint))
					.name("cutoffKeyPoints")
					.description(
						"Cutoff key points showing the interaction between and the " +
						"effection of different configuration of the cutoff arguments " +
						"`maxMixtureRate` and `minPrevalence`."
					)
				)
				.field(field -> field
					.type(GraphQLString)
					.name("assembledConsensus")
					.description(
						"Unaligned sequence consensus assembled from codon reads and " +
						"untranslated regions. Ambiguous nucleotides are included."
					)
				)
				.field(field -> field
					.type(GraphQLString)
					.name("assembledUnambiguousConsensus")
					.description(
						"Unaligned sequence consensus assembled from codon reads and " +
						"untranslated regions. Only unambiguous nucleotides are included."
					)
				)
				.field(field -> field
					.type(new GraphQLList(oBoundMutationPrevalence.get(virusName)))
					.name("mutationPrevalences")
					.argument(arg -> arg
						.type(new GraphQLList(GeneDef.enumGene.get(virusName)))
						.name("includeGenes")
						.defaultValue(Virus.getInstance(virusName).getDefaultIncludedGenes())
						.description("Genes to be included in the results")
					)
					.description("List of mutation prevalence results.")
				) 
				.field(field -> field
					.type(new GraphQLList(oAlgorithmComparison.get(virusName)))
					.name("algorithmComparison")
					.description("List of ASI comparison results.")
					.argument(aASIAlgorithmArgument.get(virusName))
					.argument(aASICustomAlgorithmArgument));
			  
			Virus<?> virus = Virus.getInstance(virusName);
			builder = virus.getVirusGraphQLExtension().extendObjectBuilder("SequenceReadsAnalysis", builder);
			return builder.build();
		}
	);

}
