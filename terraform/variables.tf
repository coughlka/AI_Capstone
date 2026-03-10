# =============================================================================
# General Variables
# =============================================================================
variable "rg" {
  description = "Resource group name."
  type        = string
  default     = "PSU-AI894-RG"
}

variable "region" {
  description = "Region of all resources"
  type        = string
  default     = "eastus"
}

# =============================================================================
# Tags
# =============================================================================
variable "tags" {
  description = "Resource tags"
  type        = map(string)
  default = {
    "Environment" = "Production"
    "Project"     = "AI894 Capstone"
  }
}

# =============================================================================
# Storage Account
# =============================================================================
variable "storage_account_name" {
  description = "Storage account name"
  type        = string
  default     = "psuai894volume"
}

# =============================================================================
# Web App
# =============================================================================
variable "app_service_plan_name" {
  description = "App service plan name"
  type        = string
  default     = "psuai894appserviceplan"
}

variable "app_service_name" {
  description = "App service name"
  type        = string
  default     = "psuai894webapp"
}

variable "acr_name" {
  description = "Azure Container Registry name"
  type        = string
  default     = "psuai894acr"
}

variable "sku_size" {
  description = "App service plan SKU size"
  type        = string
  default     = "B1"
}

# =============================================================================
# Networking
# =============================================================================
variable "vnet_address_space" {
  description = "Address space for Virtual Network"
  type        = list(string)
  default     = ["10.0.0.0/16"]
}

variable "webapp_subnet_prefix" {
  description = "Address prefix for Web App subnet"
  type        = string
  default     = "10.0.1.0/24"
}

variable "pe_subnet_prefix" {
  description = "Address prefix for Private Endpoint subnet"
  type        = string
  default     = "10.0.2.0/24"
}

